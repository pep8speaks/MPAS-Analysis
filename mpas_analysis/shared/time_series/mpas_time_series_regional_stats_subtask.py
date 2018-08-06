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

import os
import xarray
import numpy

from mpas_analysis.shared.analysis_task import AnalysisTask

from mpas_analysis.shared.io import open_mpas_dataset, write_netcdf

from mpas_analysis.shared.io.utility import build_config_full_path, \
    make_directories


class MpasTimeSeriesRegionalStatsSubtask(AnalysisTask):
    """
    Performs analysis of the time-series output of an ocean region.

    Attributes
    ----------

    mpasTimeSeriesTask : ``MpasTimeSeriesTask``
        The task that extracts the time series from MPAS monthly output

    parentTask :  ``AnalysisTask``
        The parent task, used to get the ``taskName``, ``config`` and
        ``componentName``

    timeSeriesName : str
        A name that describes the time series (e.g. a short version of
        the important field(s) in the time series) used to name the
        output file

    regionMaskSuffix : str
        The suffix of a file in the region masks folder to use to compute
        regional means.  The prefix of the file is the MPAS mesh name

    mpasVariables : list of str
        MPAS variables to include in the time series

    statsList : one or more of ['sum', 'mean']
        The stats to be computed over the region.  Sum is a cell-area-
        weighted sum.  Mean is the cell-area-weighted mean.

    refConfig : ``MpasAnalysisConfigParser``
        The configuration options for the reference run (if any)
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    def __init__(self, mpasTimeSeriesTask, parentTask, timeSeriesName,
                 regionMaskSuffix, mpasVariables, statsList=['sum', 'mean'],
                 subtaskName='timeSeriesRegionalStats', refConfig=None):
        # {{{
        """
        Construct the analysis task.

        Parameters
        ----------
        mpasTimeSeriesTask : ``MpasTimeSeriesTask``
            The task that extracts the time series from MPAS monthly output

        parentTask :  ``AnalysisTask``
            The parent task, used to get the ``taskName``, ``config`` and
            ``componentName``

        timeSeriesName : str
            A name that describes the time series (e.g. a short version of
            the important field(s) in the time series) used to name the
            output file

        regionMaskSuffix : str
            The suffix of a file in the region masks folder to use to compute
            regional means.  The prefix of the file is the MPAS mesh name

        mpasVariables : list of str
            MPAS variables to include in the time series

        statsList : one or more of ['sum', 'mean'], optional
            The stats to be computed over the region.  Sum is a cell-area-
            weighted sum.  Mean is the cell-area-weighted mean.

        subtaskName : str, optional
            The name of the subtask

        refConfig :  ``MpasAnalysisConfigParser``, optional
            Configuration options for a reference run (if any)
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        # first, call the constructor from the base class (AnalysisTask)
        super(MpasTimeSeriesRegionalStatsSubtask, self).__init__(
            config=mpasTimeSeriesTask.config,
            taskName=parentTask.taskName,
            subtaskName=subtaskName,
            componentName=parentTask.componentName,
            tags=['timeSeries'])

        self.mpasTimeSeriesTask = mpasTimeSeriesTask
        self.parentTask = parentTask
        self.timeSeriesName = timeSeriesName
        self.regionMaskSuffix = regionMaskSuffix
        self.mpasVariables = mpasVariables
        self.statsList = statsList

        parentTask.add_subtask(self)

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
        super(MpasTimeSeriesRegionalStatsSubtask, self).setup_and_check()

        self.check_analysis_enabled(
            analysisOptionName='config_am_timeseriesstatsmonthly_enable',
            raiseException=True)

        config = self.config

        mpasMeshName = config.get('input', 'mpasMeshName')
        regionMaskDirectory = config.get('regions', 'regionMaskDirectory')

        self.regionMaskFileName = '{}/{}_{}.nc'.format(
                regionMaskDirectory, mpasMeshName, self.regionMaskSuffix)

        if not os.path.exists(self.regionMaskFileName):
            raise IOError('Regional masking file {} for regional means does '
                          'not exist'.format(self.regionMaskFileName))

        # Load mesh related variables
        try:
            self.restartFileName = self.runStreams.readpath('restart')[0]
        except ValueError:
            raise IOError('No MPAS-O restart file found: need at least one '
                          'restart file for cell areas used in regional '
                          'averages.')

        # get a list of timeSeriesStats output files from the streams file,
        # reading only those that are between the start and end dates
        self.startDate = config.get('timeSeries', 'startDate')
        self.endDate = config.get('timeSeries', 'endDate')

        self.mpasTimeSeriesTask.add_variables(variableList=self.mpasVariables)

        # }}}

    def run_task(self):  # {{{
        """
        Reads time series and computes regional stats.

        Raises
        ------
        ValueError
            If an unknown stat was provided in ``self.statsList``
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        self.logger.info("\nComputing time series regional stats...")

        mpasTimeSeriesTask = self.mpasTimeSeriesTask
        config = self.config

        baseDirectory = build_config_full_path(
            config, 'output', 'timeSeriesSubdirectory')

        make_directories(baseDirectory)

        outFileName = '{}/{}.nc'.format(baseDirectory, self.timeSeriesName)

        # Load data:
        inputFile = mpasTimeSeriesTask.outputFile
        dsIn = open_mpas_dataset(fileName=inputFile,
                                 calendar=self.calendar,
                                 variableList=self.mpasVariables,
                                 startDate=self.startDate,
                                 endDate=self.endDate)
        try:
            if os.path.exists(outFileName):
                # The file already exists so load it
                dsOut = xarray.open_dataset(outFileName)
                if numpy.all(dsOut.Time.values == dsIn.Time.values):
                    return dsOut.totalMeltFlux, dsOut.meltRates
                else:
                    self.logger.warning('File {} is incomplete. Deleting '
                                        'it.'.format(outFileName))
                    os.remove(outFileName)
        except OSError:
            # something is potentailly wrong with the file, so let's delete
            # it and try again
            self.logger.warning('Problems reading file {}. Deleting '
                                'it.'.format(outFileName))
            os.remove(outFileName)

        dsRegionMask = xarray.open_dataset(self.regionMaskFileName)

        dsRegionMask = self.customize_region_masks(dsRegionMask)
        cellMasks = dsRegionMask.regionCellMasks

        dsIn, mpasVariables = self.customize_time_series_before_stats(
                dsIn, self.mpasVariables)

        dsOut = xarray.Dataset()

        dsOut['areaSum'] = (cellMasks*dsIn.areaCell).sum(dim='nCells')

        for varName in mpasVariables:
            var = dsIn[varName]
            varSum = (cellMasks*dsIn.areaCell*var).sum(dim='nCells')
            for stat in self.statsList:
                if stat == 'sum':
                    outVarName = '{}Sum'.format(varName)
                    dsOut[outVarName] = varSum
                elif stat == 'mean':
                    outVarName = '{}Mean'.format(varName)
                    dsOut[outVarName] = varSum / dsOut.areaSum
                else:
                    raise ValueError('Unknown stat {}.'.format(stat))

        dsOut = self.customize_time_series_after_stats(dsOut)

        write_netcdf(dsOut, outFileName)

        # }}}

    def customize_time_series_before_stats(self, dsTimeSeries, mpasVariables):
        # {{{
        """
        Override this function to customize the time series before the
        regional stats are computed (e.g. to apply a mask).  By default,
        adds areaCell to the data set.

        Parameters
        ----------
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

        with xarray.open_dataset(self.restartFileName) as dsRestart:
            varsToDrop = list(dsRestart.data_vars.keys())
            varsToDrop.remove('areaCell')
            dsRestart = dsRestart.drop(varsToDrop)
            dsRestart.load()
            dsTimeSeries['areaCell'] = dsRestart.areaCell

        return dsTimeSeries, mpasVariables  # }}}

    def customize_region_masks(self, dsRegionMasks):  # {{{
        """
        Override this function to customize the region masks (e.g. to
        drop some regions).

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

        return dsRegionMasks  # }}}

    def customize_time_series_after_stats(self, dsTimeSeries):  # {{{
        """
        Override this function to customize the time series after regional
        stats have been computed

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

        return dsTimeSeries  # }}}

# }}}

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
