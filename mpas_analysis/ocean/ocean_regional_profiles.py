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


import xarray
import numpy

from mpas_analysis.shared.analysis_task import AnalysisTask

from mpas_analysis.shared.constants import constants

from mpas_analysis.shared.io.utility import build_config_full_path, \
    make_directories

from mpas_analysis.shared.io import write_netcdf

from mpas_analysis.shared.climatology import compute_climatology

from mpas_analysis.shared.time_series import MpasTimeSeriesRegionalStatsSubtask


class OceanRegionalProfiles(AnalysisTask):
    # {{{
    """
    Compute and plot regionally averaged profiles of 3D MPAS fields

    Attributes
    ----------

    fields : list of dict
        dicts describing each MPAS variables, including: 'prefix', 'mpas',
        'units', 'titleName'


    seasons : list of str
        A list of seasons (keys in ``shared.constants.monthDictionary``)
        over which the climatology should be computed
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    def __init__(self, config, mpasTimeSeriesTask, refConfig=None):
        # {{{
        """
        Construct the analysis task.

        Parameters
        ----------
        config :  ``MpasAnalysisConfigParser``
            Configuration options

        mpasClimatologyTask : ``MpasClimatologyTask``
            The task that produced the climatology to be remapped and plotted

        refConfig :  ``MpasAnalysisConfigParser``, optional
            Configuration options for a reference run (if any)
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        taskName = 'oceanRegionalProfiles'
        sectionName = taskName

        self.fields = config.getExpression(sectionName, 'fields')

        prefixes = [field['prefix'] for field in self.fields]

        # call the constructor from the base class (AnalysisTask)
        super(OceanRegionalProfiles, self).__init__(
                config=config, taskName='oceanRegionalProfiles',
                componentName='ocean',
                tags=['climatology', 'profile'] + prefixes)

        config = self.config

        # read in what seasons we want to plot
        seasons = config.getExpression(sectionName, 'seasons')

        if len(seasons) == 0:
            raise ValueError('config section {} does not contain valid list '
                             'of seasons'.format(sectionName))
        self.seasons = seasons

        regionMaskSuffix = config.get(sectionName, 'regionMaskSuffix')

        self.climatologyRegionalProfileSubtask = \
            ClimatologyRegionalProfileSubtask(
                    mpasTimeSeriesTask, self,
                    timeSeriesName='regionalMeanProfiles',
                    regionMaskSuffix=regionMaskSuffix,
                    subtaskName='climatologyRegionalProfile',
                    refConfig=refConfig)

    def run_task(self):  # {{{
        """
        Plot the regional profiles
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        return
        # }}}

    # }}}


class ClimatologyRegionalProfileSubtask(MpasTimeSeriesRegionalStatsSubtask):
    # {{{
    """
    An analysis subtask computing climatologies of regional mean profiles of
    3D MPAS fields
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    def __init__(self, mpasTimeSeriesTask, parentTask, timeSeriesName,
                 regionMaskSuffix, subtaskName='climatologyRegionalProfile',
                 refConfig=None):
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

        subtaskName : str, optional
            The name of the subtask

        refConfig :  ``MpasAnalysisConfigParser``, optional
            Configuration options for a reference run (if any)
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        self.mpasFields = parentTask.fields
        self.seasons = parentTask.seasons

        mpasVariables = [field['mpas'] for field in self.mpasFields]

        # call the constructor from the base class
        # (MpasTimeSeriesRegionalStatsSubtask)
        super(ClimatologyRegionalProfileSubtask, self).__init__(
            mpasTimeSeriesTask=mpasTimeSeriesTask,
            parentTask=parentTask,
            timeSeriesName=timeSeriesName,
            regionMaskSuffix=regionMaskSuffix,
            mpasVariables=mpasVariables,
            statsList=['mean'],
            subtaskName=subtaskName,
            refConfig=refConfig)

        self.tags = self.tags + ['climatology']
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

        # first, compute regional stats by calling the base class's run_task()
        super(ClimatologyRegionalProfileSubtask, self).run_task()

        # next, compute the desired climatologies
        baseDirectory = build_config_full_path(
            self.config, 'output', 'timeSeriesSubdirectory')

        timeSeriesFileName = '{}/{}.nc'.format(
                baseDirectory, self.timeSeriesName)

        dsTimeSeries = xarray.open_dataset(timeSeriesFileName)
        for season in self.seasons:
            months = constants.monthDictionary[season]
            climatology = compute_climatology(dsTimeSeries, months,
                                              calendar=self.calendar,
                                              maskVaries=False)
            for field in self.mpasFields:
                varName = field['mpas']
                meanName = '{}Mean'.format(varName)
                stdName = '{}StdDev'.format(varName)
                climatology[stdName] = \
                    climatology['{}Squared'.format(varName)] - \
                    climatology[meanName]**2
                climatology[meanName].attrs['units'] = field['units']
                climatology[stdName].attrs['units'] = field['units']
            fileName = self.get_file_name(season)
            write_netcdf(climatology, fileName)
        # }}}

    def customize_time_series_before_stats(self, dsTimeSeries, mpasVariables):
        # {{{
        """
        Mask out values below maxLevelCell

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

        super(ClimatologyRegionalProfileSubtask,
              self).customize_time_series_before_stats(dsTimeSeries,
                                                       mpasVariables)

        with xarray.open_dataset(self.restartFileName) as dsRestart:
            nVertLevels = dsRestart.sizes['nVertLevels']
            varsToDrop = list(dsRestart.data_vars.keys())
            varsToDrop.remove('maxLevelCell')
            dsRestart = dsRestart.drop(varsToDrop)
            dsRestart.load()
            maxLevelCell = dsRestart.maxLevelCell

        vertIndex = \
            xarray.DataArray.from_dict({'dims': ('nVertLevels',),
                                        'data': numpy.arange(nVertLevels)})

        for varName in mpasVariables:
            dsTimeSeries[varName] = dsTimeSeries[varName].where(
                    vertIndex < maxLevelCell)

        return dsTimeSeries, mpasVariables  # }}}

    def customize_time_series_after_stats(self, dsTimeSeries):  # {{{
        """
        Add squared fields in addition to the fields themselves, allowing us
        to compute standard deviations in time.

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

        for varName in self.mpasVariables:
            dsTimeSeries['{}Squared'.format(varName)] = \
                dsTimeSeries['{}Mean'.format(varName)]**2

        return dsTimeSeries  # }}}

    def get_file_name(self, season):  # {{{
        """
        Get the file name for a profile climatology

        Parameters
        ----------
        season : str
            One of the seasons in ``constants.monthDictionary``
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        config = self.config
        componentName = self.componentName
        startYear = config.getint('climatology', 'startYear')
        endYear = config.getint('climatology', 'endYear')

        if componentName == 'ocean':
            ncclimoModel = 'mpaso'
        elif componentName == 'seaIce':
            ncclimoModel = 'mpascice'
        else:
            raise ValueError('component {} is not supported by ncclimo.\n'
                             'Check with Charlie Zender and Xylar Asay-Davis\n'
                             'about getting it added'.format(componentName))

        climatologyBaseDirectory = build_config_full_path(
            config, 'output', 'mpasClimatologySubdirectory')

        directory = '{}/profiles'.format(climatologyBaseDirectory)

        make_directories(directory)

        make_directories(directory)
        monthValues = sorted(constants.monthDictionary[season])
        startMonth = monthValues[0]
        endMonth = monthValues[-1]

        suffix = '{:04d}{:02d}_{:04d}{:02d}_climo'.format(
                startYear, startMonth, endYear, endMonth)

        if season in constants.abrevMonthNames:
            season = '{:02d}'.format(monthValues[0])
        fileName = '{}/{}_{}_{}.nc'.format(directory, ncclimoModel,
                                           season, suffix)
        return fileName  # }}}

    # }}}

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
