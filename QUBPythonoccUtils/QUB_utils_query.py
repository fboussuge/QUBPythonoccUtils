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

def getEntities(cursor): #get all entities having Dim>0, i.e. Body, Face and Edge
  '''[Query to get all entities having Dim>0, i.e. Body, Face and Edge]
  
  Arguments:
    cursor {[type]} -- [db cursor]
  
  Returns:
    [type] -- [tuple (entityID, )]
  '''
  r = cursor.execute(
            """
           SELECT * 
           FROM Entity
           WHERE (((Entity.Dim)>0) )
            """
                )
  ret = []
  for c in r:
    ret.append(c)
  return ret

def getFaceProperties(cursor): #get the properties of all faces: facetype
  '''[Query to get the properties of all faces: facetype]
  
  Arguments:
    cursor {[type]} -- [db cursor]
  
  Returns:
    [type] -- [tuple (FaceID, intFaceType)]
  '''
  r = cursor.execute(
            """
           SELECT Face,Type 
           FROM FaceProperties
            """
                )
  ret = []
  for c in r:
    ret.append(c)
  return ret

def getEdgeProperties(cursor): #get the properties of all edges: edgetype, edgeconvexity
  '''[Query to get the properties of all edges: edgetype, edgeconvexity]
  
  Arguments:
    cursor {[type]} -- [db cursor]
  
  Returns:
    [type] -- [tuple (EdgeID, intEdgeType, intEdgeConvexity)]
  '''
  r = cursor.execute(
            """
           SELECT * 
           FROM EdgeProperties
            """
                )
  ret = []
  for c in r:
    ret.append(c)
  return ret

def getEntitiesSolids(cursor): #get the properties of all solids (none yet)
  '''[Query to get the properties of all solids (none yet)]
  
  Arguments:
    cursor {[type]} -- [db cursor]
  
  Returns:
    [type] -- [tuple (SolidID)]
  '''
  r = cursor.execute(
            """
           SELECT * 
           FROM Entity
           WHERE (((Entity.Dim)=3) )
            """
                )
  ret = []
  for c in r:
    ret.append(c[0])
  return ret

def getBoundedFacesofBody(cursor, body): #get the bounded entities of an entity, e.g. the faces bounding a solid
  """[Query to get the bounded entities of an entity, e.g. the faces bounding a solid]
  
  Arguments:
    cursor {[type]} -- [db cursor]
    body {[type]} -- [ID of the bounding entity, e.g. solid ID]
  
  Returns:
    [type] -- [list of bounded entity ID]
  """

  sql = ''' SELECT Topology.BoundEntity
          FROM Topology INNER JOIN Entity ON Topology.Entity = Entity.Label
		  WHERE (((Entity.Label)=={0}) )'''.format(body)

  r = cursor.execute(sql)
  ret = []
  for c in r:
    ret.append(c[0])
  return ret  