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

def getEntities_VT(cursor): #get all entities which are in Topology_VT having Dim>0, i.e. Body, Face and Edge
  '''[Query to get all entities having Dim>0, i.e. Body, Face and Edge]
  
  Arguments:
    cursor {[type]} -- [db cursor]
  
  Returns:
    [type] -- [tuple (entityID, )]
  '''
  r = cursor.execute(
            """
           SELECT *
           FROM Entity INNER JOIN Topology_VT on Entity.Label = Topology_VT.Entity
           WHERE (((Entity.Dim)>0) )
           GROUP BY (Entity)
            """
                )
  ret = []
  for c in r:
    ret.append(c)
  return ret

def getFaceProperties_VT(cursor): #get the properties of all faces and superset faces: facetype 
  '''[Query to get the properties of all faces and supersets faces: facetype]
  
  Arguments:
    cursor {[type]} -- [db cursor]
  
  Returns:
    [type] -- [tuple (FaceID, intFaceType)]
  '''
  #firt query for getting original face
  r = cursor.execute(
            """
            SELECT Face, Type
            FROM FaceProperties
            WHERE Face IN (SELECT Entity FROM Topology_VT) 
            """
                )
  ret = []
  for c in r:
    ret.append(c)

  #second query for superset faces
  r = cursor.execute(
            """
            SELECT VirtualTopology.Entity, FaceProperties.Type
            FROM VirtualTopology INNER JOIN FaceProperties ON VirtualTopology.HostEntity = FaceProperties.Face
            GROUP BY VirtualTopology.Entity
            HAVING (((VirtualTopology.Entity) In (SELECT VirtualTopology.Entity
            FROM Entity INNER JOIN VirtualTopology ON Entity.Label = VirtualTopology.HostEntity
            GROUP BY Entity.Dim, VirtualTopology.Entity
            HAVING (((Entity.Dim)=2) AND ((Count(VirtualTopology.Entity))>1))))) 
            """
                )

  for c in r:
    ret.append(c)

  
  return ret

def getEdgeProperties_VT(cursor): #get the properties of all edges and superset edges: facetype
  '''[Query to get the properties of all edges and supersets edges: facetype]
  
  Arguments:
    cursor {[type]} -- [db cursor]
  
  Returns:
    [type] -- [tuple (FaceID, intFaceType)]
  '''
  #firt query for getting original edges
  r = cursor.execute(
            """
            SELECT Edge, Convexity, Type
            FROM EdgeProperties
            WHERE Edge IN (SELECT Entity FROM Topology_VT) 
            """
                )
  ret = []
  for c in r:
    ret.append(c)

  #second query for superset faces
  r = cursor.execute(
            """
            SELECT VirtualTopology.Entity, EdgeProperties.Convexity, EdgeProperties.Type
            FROM VirtualTopology INNER JOIN EdgeProperties ON VirtualTopology.HostEntity = EdgeProperties.Edge
            GROUP BY VirtualTopology.Entity
            HAVING (((VirtualTopology.Entity) In (SELECT VirtualTopology.Entity
            FROM Entity INNER JOIN VirtualTopology ON Entity.Label = VirtualTopology.HostEntity
            GROUP BY Entity.Dim, VirtualTopology.Entity
            HAVING (((Entity.Dim)=1) AND ((Count(VirtualTopology.Entity))>1))))) 
            """
                )

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
           SELECT Entity.Label, Entity.Dim 
           FROM Entity
           WHERE (((Entity.Dim)=3) )
            """
                )
  ret = []
  for c in r:
    ret.append(c)
  return ret

def getBoundedFacesofBody(cursor, body): #get the bounded faces of a solid