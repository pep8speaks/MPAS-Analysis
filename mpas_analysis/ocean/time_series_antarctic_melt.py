# This software is open source software available under the BSD-3 license.
#
# Copyright (c) 2018 Los Alamos National Security, LLC. All rights reserved.
# Copyright (c) 2018 Lawrence Livermore National Security, LLC. All rights
# reserved.
# Copyright (c) 2018 UT-Battelle, LLC. All rights reserved.
#
# Additional copyright and license information can be found in the LICENSE file
# distributed with this code, or at
# https://raw.githubusercontent.com/MPAS-Dev/MPAS-Analysis/master/LICENSE
from __future__ import absolute_import, division, print_function, \
    unicode_literals

import xarray

from mpas_analysis.shared.analysis_task import AnalysisTask

from mpas_analysis.shared.constants import constants

from mpas_analysis.shared.plot.plotting import timeseries_analysis_plot

from mpas_analysis.shared.io.utility import build_config_full_path, \
    make_directories

from mpas_analysis.shared.time_series import MpasTimeSeriesRegionalStatsSubtask

from mpas_analysis.shared.html import write_image_xml

import csv


class TimeSeriesAntarcticMelt(AnalysisTask):
    """
    Performs analysis of the time-series output of Antarctic sub-ice-shelf
    melt rates.

    Attributes
    ----------

    mpasTimeSeriesTask : ``MpasTimeSeriesTask``
        The task that extracts the time series from MPAS monthly output

    refConfig : ``MpasAnalysisConfigParser``
        The configuration options for the reference run (if any)
    """
    # Authors
    # -------
    # Xylar Asay-Davis, Stephen Price

    def __init__(self, config, mpasTimeSeriesTask, refConfig=None):
        # {{{
        """
        Construct the analysis task.

        Parameters
        ----------
        config :  ``MpasAnalysisConfigParser``
            Configuration options

        mpasTimeSeriesTask : ``MpasTimeSeriesTask``
            The task that extracts the time series from MPAS monthly output

        refConfig :  ``MpasAnalysisConfigParser``, optional
            Configuration options for a reference run (if any)
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        # first, call the constructor from the base class (AnalysisTask)
        super(TimeSeriesAntarcticMelt, self).__init__(
            config=config,
            taskName='timeSeriesAntarcticMelt',
            componentName='ocean',
            tags=['timeSeries', 'melt', 'landIceCavities'])

        self.mpasTimeSeriesTask = mpasTimeSeriesTask

        self.regionalStatsSubtask = RegionalStatsAntarcticMelt(
                mpasTimeSeriesTask, self,
                timeSeriesName='iceShelfAggregatedFluxes',
                regionMaskSuffix='iceShelfMasks',
                mpasVariables=['timeMonthly_avg_landIceFreshwaterFlux'],
                statsList=['sum', 'mean'],
                refConfig=refConfig)

        self.run_after(mpasTimeSeriesTask)
        self.refConfig = refConfig

        # }}}

    def setup_and_check(self):  # {{{
        """
        Perform steps to set up the analysis and check for errors in the setup.

        Raises
        ------
        IOError
            If files are not present
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        # first, call setup_and_check from the base class (AnalysisTask),
        # which will perform some common setup, including storing:
        #   self.inDirectory, self.plotsDirectory, self.namelist, self.streams
        #   self.calendar
        super(TimeSeriesAntarcticMelt, self).setup_and_check()

        config = self.config

        landIceFluxMode = self.namelist.get('config_land_ice_flux_mode')
        if landIceFluxMode not in ['standalone', 'coupled']:
            raise ValueError('*** timeSeriesAntarcticMelt requires '
                             'config_land_ice_flux_mode \n'
                             '    to be standalone or coupled.  Otherwise, no '
                             'melt rates are available \n'
                             '    for plotting.')

        iceShelvesToPlot = config.getExpression('timeSeriesAntarcticMelt',
                                                'iceShelvesToPlot')

        mpasMeshName = config.get('input', 'mpasMeshName')
        regionMaskDirectory = config.get('regions', 'regionMaskDirectory')

        regionMaskFileName = '{}/{}_{}.nc'.format(
                regionMaskDirectory, mpasMeshName,
                self.regionalStatsSubtask.regionMaskSuffix)

        with xarray.open_dataset(regionMaskFileName) as dsRegionMask:
            regionNames = [bytes.decode(name) for name in
                           dsRegionMask.regionNames.values]
            nRegions = dsRegionMask.dims['nRegions']

        if 'all' in iceShelvesToPlot:
            iceShelvesToPlot = regionNames
            regionIndices = [iRegion for iRegion in range(nRegions)]

        else:
            regionIndices = []
            for regionName in iceShelvesToPlot:
                if regionName not in regionNames:
                    raise ValueError('Unknown ice shelf name {}'.format(
                            regionName))

                iRegion = regionNames.index(regionName)
                regionIndices.append(iRegion)

        self.regionIndices = regionIndices
        self.iceShelvesToPlot = iceShelvesToPlot
        self.xmlFileNames = []

        for prefix in ['melt_flux', 'melt_rate']:
            for regionName in iceShelvesToPlot:
                regionName = regionName.replace(' ', '_')
                self.xmlFileNames.append(
                    '{}/{}_{}.xml'.format(self.plotsDirectory, prefix,
                                          regionName))
        return  # }}}

    def run_task(self):  # {{{
        """
        Performs analysis of the time-series output of Antarctic sub-ice-shelf
        melt rates.
        """
        # Authors
        # -------
        # Xylar Asay-Davis, Stephen Price

        self.logger.info("\nPlotting Antarctic melt rate time series...")

        self.logger.info('  Load melt rate data...')

        config = self.config
        calendar = self.calendar

        totalMeltFlux, meltRates = \
            self._load_ice_shelf_fluxes(self.config)

        plotRef = self.refConfig is not None
        if plotRef:
            refRunName = self.refConfig.get('runs', 'mainRunName')

            refTotalMeltFlux, refMeltRates = \
                self._load_ice_shelf_fluxes(self.refConfig)

        # Load observations from multiple files and put in dictionary based
        # on shelf keyname
        observationsDirectory = build_config_full_path(config,
                                                       'oceanObservations',
                                                       'meltSubdirectory')
        obsFileNameDict = {'Rignot et al. (2013)':
                           'Rignot_2013_melt_rates.csv',
                           'Rignot et al. (2013) SS':
                           'Rignot_2013_melt_rates_SS.csv'}

        obsDict = {}  # dict for storing dict of obs data
        for obsName in obsFileNameDict:
            obsFileName = '{}/{}'.format(observationsDirectory,
                                         obsFileNameDict[obsName])
            obsDict[obsName] = {}
            obsFile = csv.reader(open(obsFileName, 'rU'))
            next(obsFile, None)  # skip the header line
            for line in obsFile:  # some later useful values commented out
                shelfName = line[0]
                # surveyArea = line[1]
                meltFlux = float(line[2])
                meltFluxUncertainty = float(line[3])
                meltRate = float(line[4])
                meltRateUncertainty = float(line[5])
                # actualArea = float( line[6] )  # actual area here is in sq km

                # build dict of obs. keyed to filename description
                # (which will be used for plotting)
                obsDict[obsName][shelfName] = {
                        'meltFlux': meltFlux,
                        'meltFluxUncertainty': meltFluxUncertainty,
                        'meltRate': meltRate,
                        'meltRateUncertainty': meltRateUncertainty}

        # If areas from obs file used need to be converted from sq km to sq m

        mainRunName = config.get('runs', 'mainRunName')
        movingAverageMonths = config.getint('timeSeriesAntarcticMelt',
                                            'movingAverageMonths')

        nRegions = totalMeltFlux.sizes['nRegions']

        outputDirectory = build_config_full_path(config, 'output',
                                                 'timeseriesSubdirectory')

        make_directories(outputDirectory)

        self.logger.info('  Make plots...')
        for iRegion in range(nRegions):

            regionName = self.iceShelvesToPlot[iRegion]

            # get obs melt flux and unc. for shelf (similar for rates)
            obsMeltFlux = []
            obsMeltFluxUnc = []
            obsMeltRate = []
            obsMeltRateUnc = []
            for obsName in obsDict:
                if regionName in obsDict[obsName]:
                    obsMeltFlux.append(
                        obsDict[obsName][regionName]['meltFlux'])
                    obsMeltFluxUnc.append(
                        obsDict[obsName][regionName]['meltFluxUncertainty'])
                    obsMeltRate.append(
                        obsDict[obsName][regionName]['meltRate'])
                    obsMeltRateUnc.append(
                        obsDict[obsName][regionName]['meltRateUncertainty'])
                else:
                    # append NaN so this particular obs won't plot
                    self.logger.warning('{} observations not available for '
                                        '{}'.format(obsName, regionName))
                    obsMeltFlux.append(None)
                    obsMeltFluxUnc.append(None)
                    obsMeltRate.append(None)
                    obsMeltRateUnc.append(None)

            title = regionName.replace('_', ' ')

            regionName = regionName.replace(' ', '_')

            xLabel = 'Time (yr)'
            yLabel = 'Melt Flux (GT/yr)'

            timeSeries = totalMeltFlux.isel(nRegions=iRegion)

            filePrefix = 'melt_flux_{}'.format(regionName)
            figureName = '{}/{}.png'.format(self.plotsDirectory, filePrefix)

            fields = [timeSeries]
            lineColors = ['k']
            lineWidths = [2.5]
            legendText = [mainRunName]
            if plotRef:
                fields.append(refTotalMeltFlux.isel(nRegions=iRegion))
                lineColors.append('r')
                lineWidths.append(1.2)
                legendText.append(refRunName)

            timeseries_analysis_plot(config, fields, movingAverageMonths,
                                     title, xLabel, yLabel, figureName,
                                     calendar=calendar,
                                     lineColors=lineColors,
                                     lineWidths=lineWidths,
                                     legendText=legendText,
                                     obsMean=obsMeltFlux,
                                     obsUncertainty=obsMeltFluxUnc,
                                     obsLegend=list(obsDict.keys()))

            caption = 'Running Mean of Total Melt Flux  under Ice ' \
                      'Shelves in the {} Region'.format(title)
            write_image_xml(
                config=config,
                filePrefix=filePrefix,
                componentName='Ocean',
                componentSubdirectory='ocean',
                galleryGroup='Antarctic Melt Time Series',
                groupLink='antmelttime',
                gallery='Total Melt Flux',
                thumbnailDescription=title,
                imageDescription=caption,
                imageCaption=caption)

            xLabel = 'Time (yr)'
            yLabel = 'Melt Rate (m/yr)'

            timeSeries = meltRates.isel(nRegions=iRegion)

            filePrefix = 'melt_rate_{}'.format(regionName)
            figureName = '{}/{}.png'.format(self.plotsDirectory, filePrefix)

            fields = [timeSeries]
            lineColors = ['k']
            lineWidths = [2.5]
            legendText = [mainRunName]
            if plotRef:
                fields.append(refMeltRates.isel(nRegions=iRegion))
                lineColors.append('r')
                lineWidths.append(1.2)
                legendText.append(refRunName)

            if config.has_option(self.taskName, 'firstYearXTicks'):
                firstYearXTicks = config.getint(self.taskName,
                                                'firstYearXTicks')
            else:
                firstYearXTicks = None

            if config.has_option(self.taskName, 'yearStrideXTicks'):
                yearStrideXTicks = config.getint(self.taskName,
                                                 'yearStrideXTicks')
            else:
                yearStrideXTicks = None

            timeseries_analysis_plot(config, fields, movingAverageMonths,
                                     title, xLabel, yLabel, figureName,
                                     calendar=calendar,
                                     lineColors=lineColors,
                                     lineWidths=lineWidths,
                                     legendText=legendText,
                                     obsMean=obsMeltRate,
                                     obsUncertainty=obsMeltRateUnc,
                                     obsLegend=list(obsDict.keys()),
                                     firstYearXTicks=firstYearXTicks,
                                     yearStrideXTicks=yearStrideXTicks)

            caption = 'Running Mean of Area-averaged Melt Rate under Ice ' \
                      'Shelves in the {} Region'.format(title)
            write_image_xml(
                config=config,
                filePrefix=filePrefix,
                componentName='Ocean',
                componentSubdirectory='ocean',
                galleryGroup='Antarctic Melt Time Series',
                groupLink='antmelttime',
                gallery='Area-averaged Melt Rate',
                thumbnailDescription=title,
                imageDescription=caption,
                imageCaption=caption)
        # }}}

    def _load_ice_shelf_fluxes(self, config):  # {{{
        """
        Reads melt flux time series and computes regional total melt flux and
        mean melt rate.
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        baseDirectory = build_config_full_path(
            config, 'output', 'timeSeriesSubdirectory')

        outFileName = '{}/{}.nc'.format(
                baseDirectory, self.regionalStatsSubtask.timeSeriesName)

        dsOut = xarray.open_dataset(outFileName)
        return dsOut.totalMeltFlux, dsOut.meltRates  # }}}

# }}}


class RegionalStatsAntarcticMelt(MpasTimeSeriesRegionalStatsSubtask):  # {{{
    """
    A class for customizing the regional stats for Antarctic melt rates and
    melt fluxes.
    """

    def customize_time_series_before_stats(self, dsTimeSeries, mpasVariables):
        # {{{
        """
        Customize the time series by masking ``areaCell`` by
        ``landIceFraction``

        Parameters
        dsTimeSeries
        dsTimeSeries : ``xarray.Dataset```
            The MPAS time series data set before regional stats have been
            computed

        mpasVariables : list of str
            MPAS variables that were included in the time series

        Returns
        -------
        dsTimeSeries : ``xarray.Dataset```
            The same data set with any custom fields added or modifications
            made

        mpasVariables : list of str
            A possibly modified list of MPAS variables on which to compute
            stats

        """
        # Authors
        # -------
        # Xylar Asay-Davis

        super(RegionalStatsAntarcticMelt,
              self).customize_time_series_before_stats(dsTimeSeries,
                                                       mpasVariables)

        with xarray.open_dataset(self.restartFileName) as dsRestart:
            varsToDrop = list(dsRestart.data_vars.keys())
            varsToDrop.remove('landIceMask')
            dsRestart = dsRestart.drop(varsToDrop)
            dsRestart.load()
            landIceFraction = dsRestart.landIceFraction.isel(Time=0)

            dsTimeSeries['areaCell'] = landIceFraction*dsTimeSeries.areaCell

        return dsTimeSeries, mpasVariables  # }}}

    def customize_region_masks(self, dsRegionMasks):  # {{{
        """
        Cusomize the region mask to use only the requested ice shelf regions
        rather than all ice shelves

        Parameters
        ----------
        dsRegionMasks : ``xarray.Dataset```
            The MPAS region masks data set

        Returns
        -------
        dsRegionMasks : ``xarray.Dataset```
            The same data set with any custom fields added or modifications
            made
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        dsRegionMasks = dsRegionMasks.isel(
                nRegions=self.parentTask.regionIndices)

        return dsRegionMasks  # }}}

    def customize_time_series_after_stats(self, dsTimeSeries):  # {{{
        """
        Modify the units of the resulting melt fluxes and melt rates.  Rename
        the fields to more user-friendly names.  Add units and descriptions
        for the fields.

        Parameters
        ----------
        dsTimeSeries : ``xarray.Dataset```
            The MPAS time series data set after regional stats have been
            computed

       Returns
        -------
        dsTimeSeries : ``xarray.Dataset```
            The same data set with any custom fields added or modifications
            made
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        dsTimeSeries = dsTimeSeries.rename(
                {'timeMonthly_avg_landIceFreshwaterFluxSum': 'totalMeltFlux',
                 'timeMonthly_avg_landIceFreshwaterFluxMean': 'meltRates'})

        # convert from kg/s to GT/yr
        dsTimeSeries['totalMeltFlux'] = \
            constants.sec_per_year / constants.kg_per_GT * \
            dsTimeSeries.totalMeltFlux

        # from kg/m^2/yr to m/yr
        dsTimeSeries['meltRates'] = \
            constants.sec_per_year / constants.rho_fw * \
            dsTimeSeries.meltRates

        dsTimeSeries.totalMeltFlux.attrs['units'] = 'GT a$^{-1}$'
        dsTimeSeries.totalMeltFlux.attrs['description'] = \
            'Total melt flux summed over each ice shelf or region'
        dsTimeSeries.meltRates.attrs['units'] = 'm a$^{-1}$'
        dsTimeSeries.meltRates.attrs['description'] = \
            'Melt rate averaged over each ice shelf or region'

        return dsTimeSeries  # }}}
    # }}}

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
