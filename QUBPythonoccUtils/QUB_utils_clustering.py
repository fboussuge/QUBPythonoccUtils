# coding: utf-8
# Copyright (C) 2018 Flavien Boussuge <f.boussuge@qub.ac.uk>, Liang Sun <Liang.Sun@qub.ac.uk>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""
QUBPythonoccUtils Python package

This package provides functions/methods:
  - to query to an sql database containing data on CAD models

"""
from matplotlib import pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.cluster.hierarchy import fcluster

def full_dendrogram(*args, **kwargs): #create a full dendogram graph
    """[create a full dendogram graph]
    
    Returns:
        [type] -- [dendogram to plot]
    """
    #calculate the full dendrogram
    ddata = dendrogram(*args, **kwargs)
    plt.title('Hierarchical Clustering Dendogram')
    plt.xlabel('sample index')
    plt.ylabel('distance')
    #plt.yscale('log')
    #plt.ylim([1e-20, 1.0])
    return ddata

def custom_dendrogram(*args, **kwargs): #create a custom dendogram graph
    """[create a custom dendogram graph]
    
    Returns:
        [type] -- [dendogram to plot]
    """

    max_d = kwargs.pop('max_d', None)
    if max_d and 'color_threshold' not in kwargs:
        kwargs['color_threshold'] = max_d
    annotate_above = kwargs.pop('annotate_above', 0)
    
    ddata = dendrogram(*args, **kwargs)
    
    if not kwargs.get('no_plot', False):
        plt.title('Hierarchical Clustering Dendrogram (truncated)')
        plt.xlabel('sample index or (cluster size)')
        plt.ylabel('distance')
        # plt.yscale('log')
        # plt.ylim([1e-20, 1.0])
        for i, d, c in zip(ddata['icoord'], ddata['dcoord'], ddata['color_list']):
            x = 0.5 * sum(i[1:3])
            y = d[1]
            if y > annotate_above:
                plt.plot(x, y, 'o', c=c)
                plt.annotate("%.3g" % y, (x,y), xytext=(0, -5),
                             textcoords='offset points',
                             va='top', ha='center')
        if max_d:
            plt.axhline(y=max_d, c='k')
    return ddata
