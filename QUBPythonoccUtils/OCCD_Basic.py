# -*- coding: utf-8 -*-
"""
Created on Tue Jan  9 15:07:47 2018

@author: 3052180 Liang Sun
useful website for OCCD:
====================================
file:///C:/pythonocc-core-0.18/doc/apidoc/0.18/search.html?q=D1&check_keywords=yes&area=default
https://www.opencascade.com/doc/occt-6.9.1/refman/html/class_topo_d_s___shape.html
https://www.opencascade.com/content/point-parameter
====================================
"""
import sys
import numpy as np

from OCC.gp import gp_Pnt, gp_Vec
from OCC.ShapeAnalysis import ShapeAnalysis_Surface
from OCC.BRep import BRep_Tool
from OCC.STEPControl import STEPControl_Reader
from OCC.IFSelect import IFSelect_RetDone, IFSelect_ItemsByEntity
from OCC.GeomLProp import GeomLProp_SLProps 
from OCC.TopAbs import TopAbs_FORWARD,TopAbs_REVERSED
from OCC.BRepTools import breptools
from OCC.BRepAdaptor import BRepAdaptor_Curve
from OCC.GCPnts import GCPnts_AbscissaPoint, GCPnts_AbscissaPoint_Length
from OCC.BRepBuilderAPI import BRepBuilderAPI_MakeVertex
from OCC.BRepExtrema import BRepExtrema_DistShapeShape


from OCC.TDocStd import Handle_TDocStd_Document
from OCC.XCAFApp import XCAFApp_Application
from OCC.TCollection import TCollection_ExtendedString
from OCC.XCAFDoc import (XCAFDoc_DocumentTool_ShapeTool,
                         XCAFDoc_DocumentTool_ColorTool,
                         XCAFDoc_DocumentTool_LayerTool,
                         XCAFDoc_DocumentTool_MaterialTool,
                         XCAFDoc_ColorSurf,XCAFDoc_ColorGen)
from OCC.STEPCAFControl import STEPCAFControl_Reader
from OCC.TDF import TDF_LabelSequence
from OCC.Quantity import Quantity_Color,Quantity_TOC_RGB
from OCC.ShapeFix import ShapeFix_Wireframe,ShapeFix_Shape
from OCCD_Topo_Traverse import Topo
from OCC.BRepGProp import brepgprop
from OCC.GProp import GProp_GProps

import Python_Basic as PB


#############################################################
####                      Topology                      #####
#############################################################



def ask_face_wires(face, topo):
    """
    This returns a list of wires(which is a loop of oriented edge). The first one is the boundary of the face.
    The followings (if has) are the inner loops.
    """
    wireList=[]
    allWires=list(topo.wires_from_face(face))
    
    if len(allWires)==1:
        return allWires
    elif len(allWires)>1:
         outer=breptools.OuterWire(face)
         inners=[x for x in allWires if x!=outer]
         wireList.append(outer)
         wireList.extend(inners)
         return wireList
    else:
        print("Function: face_wires, the number of wires of the face is 0 ")    
        return None



def get_edges_bounded_by_vertex_in_face(topo,vertex, face):
    """
    """
    boundingEdge=[]
    es1=list(topo.edges_from_vertex(vertex))
    es2=list(topo.edges_from_face(face))
    boundingEdge=PB.list_common(es1,es2)
    if len(boundingEdge)!=2:
        print("The vertex is not bounded by two edges")
    return boundingEdge


def get_wire_of_edge_in_face(topo,edge, face):
    """
    """
    ws1=list(topo.wires_from_edge(edge))
    ws2=list(topo.wires_from_face(face))
    commonW=PB.list_common(ws1,ws2)
    return commonW[0]


def order_two_edges_in_wire(topo,wire,e1,e2):
    """
    """
    orderedE=list(topo.ordered_edges_from_wire(wire))
    
    index1=-1
    index2=-1
    firstE=None
    nextE=None
    
    for i in range(len(orderedE)):
        if orderedE[i]==e1:
            index1=i
            break
    for i in range(len(orderedE)):
        if orderedE[i]==e2:
            index2=i
            break
    
    if index1==0 and index2==len(orderedE)-1:
        firstE=e2
        nextE=e1
    elif index2==0 and index1==len(orderedE)-1:
        firstE=e1
        nextE=e2
    elif index1<index2:
        firstE=e1
        nextE=e2
    else:
        firstE=e2
        nextE=e1
    return [firstE,nextE]


#############################################################
####                      Geometry                      #####
#############################################################

def ask_vertice_parm_edge(vertex, edge):
    """
    Ask the parameter of an edge
    """
    t=BRep_Tool().Parameter(vertex,edge)
    return t


def ask_vertex_uv_face(vertex, face):
    """
    This return the uv of a vertex on a face
    """
    gpPnt2D=BRep_Tool().Parameters(vertex,face)
    return gpPnt2D.Coord()



def ask_vertex_coor(vertex):
    """
    Ask the coordinate of the vertex
    """
    coor=BRep_Tool().Pnt(vertex).Coord()
    return coor



def ask_edge_midPnt(edge):
    """
    """
    result=BRep_Tool.Curve(edge)###result[0] is the handle of curve;result[1] is the umin; result[2] is umax
    tmid=(result[1]+result[2])/2
    p=gp_Pnt(0,0,0)
    result[0].GetObject().D0(tmid,p)
    return p.Coord()


def ask_edge_midpnt_tangent(edge):
    """
    Ask the midpoint of an edge and the tangent at the midpoint
    """
    result=BRep_Tool.Curve(edge)###result[0] is the handle of curve;result[1] is the umin; result[2] is umax
    tmid=(result[1]+result[2])/2
    p=gp_Pnt(0,0,0)
    v1=gp_Vec(0,0,0)
    result[0].GetObject().D1(tmid,p,v1)###handle.GetObject() gives Geom_Curve type, p:gp_Pnt, v1:gp_Vec
    return[p.Coord(),v1.Coord()]

def ask_edge_tangent(edge,parm):
    """
    parm is normalized parameter
    """
    result=BRep_Tool.Curve(edge)###result[0] is the handle of curve;result[1] is the umin; result[2] is umax
    t=result[1]+(result[2]-result[1])*parm
    p=gp_Pnt(0,0,0)
    v1=gp_Vec(0,0,0)
    result[0].GetObject().D1(t,p,v1)###handle.GetObject() gives Geom_Curve type, p:gp_Pnt, v1:gp_Vec
    return v1.Coord()


def ask_edge_tangent2(edge,unNormParm):
    """
    """
    result=BRep_Tool.Curve(edge)###result[0] is the handle of curve;result[1] is the umin; result[2] is umax

    p=gp_Pnt(0,0,0)
    v1=gp_Vec(0,0,0)
    result[0].GetObject().D1(unNormParm,p,v1)###handle.GetObject() gives Geom_Curve type, p:gp_Pnt, v1:gp_Vec
    return v1.Coord()


def ask_face_centroid(face):
    """
    """
    massProps=GProp_GProps()
    brepgprop.SurfaceProperties(face,massProps)
    gPt=massProps.CentreOfMass()
    return gPt.Coord()



def ask_point_uv(xyz, face, uvrange):
    """
    This is a general function which gives the uv coordiates from the xyz coordinates.
    The uv value is not normlized.
    """
    gpPnt=gp_Pnt(float(xyz[0]), float(xyz[1]), float(xyz[2]))
    surface=BRep_Tool().Surface(face) ### Handle_Geom_Surface

    sas=ShapeAnalysis_Surface(surface)
    gpPnt2D=sas.ValueOfUV(gpPnt,0.01)
    uv=list(gpPnt2D.Coord())
#    print(uvrange)
#    print(uv)
    geom_surface=surface.GetObject()
    tol=0.001
    if geom_surface.IsUPeriodic()==True:
#        print (geom_surface.UPeriod())
        if uv[0]<uvrange[0] and uvrange[0]-uv[0]>tol:
            uv[0]=uv[0]+geom_surface.UPeriod()
        if uv[0]>uvrange[1] and uv[0]-uvrange[1]>tol:
            uv[0]=uv[0]-geom_surface.UPeriod()

    
    if geom_surface.IsVPeriodic()==True:
        if uv[1]<uvrange[2] and uvrange[2]-uv[1]>tol:
            uv[1]=uv[1]+geom_surface.VPeriod()
        if uv[1]>uvrange[3] and uv[1]-uvrange[3]>tol:
            uv[1]=uv[1]-geom_surface.VPeriod()
#    print(uv)
    return uv





def ask_point_uv2(xyz, face):
    """
    This is a general function which gives the uv coordiates from the xyz coordinates.
    The uv value is not normlized.
    """
    gpPnt=gp_Pnt(float(xyz[0]), float(xyz[1]), float(xyz[2]))
    surface=BRep_Tool().Surface(face) ### Handle_Geom_Surface

    sas=ShapeAnalysis_Surface(surface)
    gpPnt2D=sas.ValueOfUV(gpPnt,0.01)
    uv=list(gpPnt2D.Coord())

    return uv






def ask_point_normal_face(uv,face):
    """
    Ask the normal vector of a point given the uv coordinate of the point on a face
    """
    surface=BRep_Tool().Surface(face)
    props=GeomLProp_SLProps(surface,uv[0],uv[1],1,1e-6)
#    GeomLProp_SLProps.SetParameters(surface,uv[0],uv[1])
#    GeomLProp_SLProps.SetSurface(surface)

    gpDir=props.Normal()#gp_Dir type
    if face.Orientation()==TopAbs_REVERSED:
        gpDir.Reverse()
        #print("face reversed")
    return gpDir.Coord()



def corner_angle(topo, vertex, face):
    """
    """
    angle=0
    edges=get_edges_bounded_by_vertex_in_face(topo,vertex, face)
    if edges==[]: ### it is possible that a vertex of a superface is only a vertex of
    ### one of the hosts. In this case, for the other host faces, the angle should be 0
        return angle
    wire=get_wire_of_edge_in_face(topo,edges[0], face)
    [firstE,nextE]=order_two_edges_in_wire(topo,wire,edges[0],edges[1])
    
    result1=edge_extreme(firstE)
    result2=edge_extreme(nextE)
    
    parm1=ask_vertice_parm_edge(vertex,firstE)
    parm2=ask_vertice_parm_edge(vertex,nextE)
#    t1=ask_edge_tangent2(firstE,parm1)
#    t2=ask_edge_tangent2(nextE,parm2)
#    print(result1[1:])
#    print(result2[1:])
#    print(parm1)
#    print(parm2)
#    print(t1)
#    print(t2)

    if abs(parm1-result1[1])<abs(parm1-result1[2]):
        t1=ask_edge_tangent2(firstE,result1[1])
    else:
        t1=ask_edge_tangent2(firstE,result1[2])
        t1 = [x * -1 for x in t1] 
    
    if abs(parm2-result2[1])<abs(parm2-result2[2]):
        t2=ask_edge_tangent2(nextE,result2[1])
    else:
        t2=ask_edge_tangent2(nextE,result2[2])
        t2 = [x * -1 for x in t2]
#    print(t1)
#    print(t2)
    xyz=ask_vertex_coor(vertex)
    uv=ask_point_uv(xyz, face,face_extreme(face))
    n=ask_point_normal_face(uv,face)
    norm0=np.sqrt(np.dot(t1,t1))
    norm1=np.sqrt(np.dot(t2,t2))
    dp=np.dot(t1,t2)
    angle=np.arccos(round(dp/(norm0*norm1),4))
#    print(angle)
    cp=np.cross(t2,t1)
    r=np.dot(cp,n)
    s=np.sign(r)
    if s==-1:
        angle=np.pi*2-angle
    else:
        pass
#    print(angle)
    return angle



def arc_length_t_edge(edge, t1, t2):
    """
    This calculates the arc length between two paramters t1 and t2
    """
    adaptor3d_curve=BRepAdaptor_Curve (edge)
    arclength=GCPnts_AbscissaPoint_Length(adaptor3d_curve,t1,t2)
    return arclength








def edge_length(edge):
#    result=BRep_Tool.Curve(edge)###result[0] is the handle of curve;result[1] is the umin; result[2] is umax
#    """
    """
    Calculate the lengh of an edge
    """
    length= GCPnts_AbscissaPoint().Length (BRepAdaptor_Curve (edge))
    return length

#def curve_length(curve):##Geom_Curve
#    curve.GetHandle()
#    GCPoints_AbscissaPoint.Length


def edge_extreme(edge):
    """
    Get the parameter range of a curve underlying an edge
    """
    (curve,tmin,tmax)=BRep_Tool.Curve(edge)
    return (curve,tmin,tmax)

def edge_dihedral(e, topo):
    """
    Calculate the dihedral angle of an edge
    """
    
    fs=list(topo.faces_from_edge(e))
    #print(len(fs))
    if len(fs)==1:
        return 0

    edges=topo.edges_from_face(fs[0])
    tem=None
    #orientation=0
    for ee in edges:
        if ee.IsEqual(e) or ee.IsSame(e):
            tem=ee
            break

    

    [midPnt,tangent]=ask_edge_midpnt_tangent(tem)
    uv0=ask_point_uv2(midPnt,fs[0])
    uv1=ask_point_uv2(midPnt,fs[1])
    n0=ask_point_normal_face(uv0,fs[0])
#    p0=xyz_from_uv_face(uv0,fs[0])
    n1=ask_point_normal_face(uv1,fs[1])
#    p1=xyz_from_uv_face(uv1,fs[1])


    
    angle=0
    if tem.Orientation()==TopAbs_FORWARD:
        #print("forward")
        dp=np.dot(n0,n1)
        norm0=np.sqrt(np.dot(n0,n0))
        norm1=np.sqrt(np.dot(n1,n1))
        angle=np.arccos(round(dp/norm0*norm1,4))
        cp=np.cross(n0,n1)
        r=np.dot(cp,tangent)
        s=np.sign(r)
        if s==-1:
            angle=np.pi+angle
        elif s==1:
            angle=np.pi-angle
        elif s==0:
            angle=np.pi
    else:
        #print("reversed")
        dp=np.dot(n0,n1)
        norm0=np.sqrt(np.dot(n0,n0))
        norm1=np.sqrt(np.dot(n1,n1))
        angle=np.arccos(round(dp/norm0*norm1,4))
        cp=np.cross(n1,n0)
        r=np.dot(cp,tangent)
        s=np.sign(r)
        if s==-1:
            angle=np.pi+angle
        elif s==1:
            angle=np.pi-angle
        elif s==0:
            angle=np.pi
    
#        print(angle)
#        print("=============")
    return angle


def face_extreme(face):
    """
    """
    uv=breptools.UVBounds(face)
    return uv


def min_distance(pt,face):
    """
    Minimum distance between point and face
    """

    gpPnt=gp_Pnt(pt[0], pt[1],pt[2])
    vertex = BRepBuilderAPI_MakeVertex(gpPnt).Vertex()
    dis=BRepExtrema_DistShapeShape(vertex,face)

    p=dis.PointOnShape2(1)
    d=dis.Value()
    return [d,p]


def make_vertex(pt):
    """
    """
    gpPnt=gp_Pnt(pt[0], pt[1],pt[2])
    vertex = BRepBuilderAPI_MakeVertex(gpPnt).Vertex()

    return vertex




def reparameterization_arclength(edge,t1,tmin):
    """
    This function reparameterize the curve using arc length as parameter.
    """
#    (curve,tmin,tmax)=BRep_Tool.Curve(edge)
    partLength=arc_length_t_edge(edge, t1, tmin)
    wholeLength= GCPnts_AbscissaPoint().Length(BRepAdaptor_Curve (edge))
    return partLength/wholeLength




def xyz_from_t_edge(edge,parameter):
    """
    Given an edge and a normalized parameter[0,1], this calculates the point on the edge
    1. First extract the curve from edge and non-normalized parameter range tmin, tmax
    2. Tranlate the input parameter to non-normalized parameter
    3. Calculate the point (return a gpPnt)
    
    """
   
    result=BRep_Tool.Curve(edge)###result[0] is the handle of curve;result[1] is the umin; result[2] is umax
#    print(parameter)
#    print(result[2])
#    print(result[1])
    np=result[1]+(result[2]-result[1])*parameter   
    aPnt=result[0].GetObject().Value(np)###handle_Geom_Curve.GetObject() will give the Geom_Curve
    return aPnt



def xyz_from_arclength_edge(edge,arclength, tmin):
    """
    Given arclength from the endpoint with parameter tmin, return the xyz on edge
    """
    adaptor3d_Curve=BRepAdaptor_Curve(edge)
    abscissaPoint= GCPnts_AbscissaPoint(adaptor3d_Curve,arclength, tmin)

    result=BRep_Tool.Curve(edge)###result[0] is the handle of curve;result[1] is the umin; result[2] is umax
    aPnt=result[0].GetObject().Value(abscissaPoint.Parameter())

    return aPnt


def xyz_from_uv_face(uv, face):
    """
    Given a face and a uv coordinate, return the xyz coordinate
    """
#    print(uv)
    [umin,umax,vmin,vmax]=face_extreme(face)
    u=umin+(umax-umin)*uv[0]
    v=vmin+(vmax-vmin)*uv[1]
    result=BRep_Tool.Surface(face)###result[0] is the handle of curve;result[1] is the umin; result[2] is umax 
    aPnt=result.GetObject().Value(u,v)###handle_Geom_Curve.GetObject() will give the Geom_Curve
    
    return aPnt


def xyz_from_uv_face_unnormlized(uv, face):
    """
    Given a face and a uv coordinate, return the xyz coordinate
    """
#    print(uv)
    result=BRep_Tool.Surface(face)###result[0] is the handle of curve;result[1] is the umin; result[2] is umax 
    aPnt=result.GetObject().Value(uv[0],uv[1])###handle_Geom_Curve.GetObject() will give the Geom_Curve
    
    return aPnt





#############################################################
####                     Attribute                      #####
#############################################################

def get_colour_from_step_file(face,colour_tool):
    """
    This is to get the colour of an topo shape
    """
    aColor=Quantity_Color()
    colour_tool.GetColor(face,XCAFDoc_ColorSurf,aColor)
    ### Time 255 since in the step the colour value is normizad betwen [0,1]
#    for x in aColor.Values(Quantity_TOC_RGB):
#        print (x)
#        print(str(255*x))
#        print(str(round(255*x)))
    rgb=tuple(round (255*x) for x in aColor.Values(Quantity_TOC_RGB))
    #print(rgb)
    return rgb


def read_step_file(filename):
    """ 
    Read the STEP file and returns a shape. This cann't access any colour/name/layer atribute
    If these are used, please use the following method instead
    """
    step_reader = STEPControl_Reader()
    status = step_reader.ReadFile(filename)

    if status == IFSelect_RetDone:  # check status
        failsonly = False
        step_reader.PrintCheckLoad(failsonly, IFSelect_ItemsByEntity)
        step_reader.PrintCheckTransfer(failsonly, IFSelect_ItemsByEntity)    
        step_reader.TransferRoot(1)
        step_reader.NbShapes()
        aResShape = step_reader.Shape(1)
        topo=Topo(aResShape)
    else:
        print("Error: can't read file.")
        sys.exit(0)
    return topo


def read_step_file_with_attribute(filename):
    """
    Read a step file and returns a shape (geometrical information). This uses the XDE of OCCD
    and will also be able to extract the attribute information such as colour or name or layer.
    """
    h_doc = Handle_TDocStd_Document()
    #print "Empty Doc?", h_doc.IsNull()
        
    ### Create the application
    app = XCAFApp_Application.GetApplication().GetObject()
    app.NewDocument(TCollection_ExtendedString("MDTV-CAF"),h_doc)
    
    ### Get tools, here only enable the colour attribute
    doc = h_doc.GetObject()
    h_shape_tool = XCAFDoc_DocumentTool_ShapeTool(doc.Main())
    h_colour_tool = XCAFDoc_DocumentTool_ColorTool(doc.Main())
#    l_Layers = XCAFDoc_DocumentTool_LayerTool(doc.Main())
#    l_materials = XCAFDoc_DocumentTool_MaterialTool(doc.Main())
    shape_tool = h_shape_tool.GetObject()
    colour_tool = h_colour_tool.GetObject()

    ### Read files
    STEPReader = STEPCAFControl_Reader()
    STEPReader.SetColorMode(True)
#    STEPReader.SetLayerMode(True)
#    STEPReader.SetNameMode(True)
#    STEPReader.SetMatMode(True)
    
    status = STEPReader.ReadFile(filename)
    STEPReader.Reader().PrintCheckLoad(True,0)
    if status == IFSelect_RetDone:
        STEPReader.Transfer(doc.GetHandle())
    
    ### Get root assembly
    shapeLabels = TDF_LabelSequence()
    shape_tool.GetFreeShapes(shapeLabels)
    print ('Number of shapes at root :%i'%shapeLabels.Length())
    
    ### Here we only have one solid in the model so we directly use the first solid
    ### If there are more, please use the for loop to get all solids
    shape=shape_tool.GetShape(shapeLabels.Value(1))
    
    result=fix_shape(shape)
#    for i in range(shapeLabels.Length()):
#        #print(Labels.Value(i+1).Tag())
#        shape=shape_tool.GetShape(shapeLabels.Value(i+1))
#        shapes.append(shape)
    
    return[result,colour_tool]



def read_step_file_with_attribute2(filename):
    """
    Read a step file and returns a shape (geometrical information). This uses the XDE of OCCD
    and will also be able to extract the attribute information such as colour or name or layer.
    """
    h_doc = Handle_TDocStd_Document()
    #print "Empty Doc?", h_doc.IsNull()
        
    ### Create the application
    app = XCAFApp_Application.GetApplication().GetObject()
    app.NewDocument(TCollection_ExtendedString("MDTV-CAF"),h_doc)
    
    ### Get tools, here only enable the colour attribute
    doc = h_doc.GetObject()
    h_shape_tool = XCAFDoc_DocumentTool_ShapeTool(doc.Main())
    h_colour_tool = XCAFDoc_DocumentTool_ColorTool(doc.Main())
#    l_Layers = XCAFDoc_DocumentTool_LayerTool(doc.Main())
#    l_materials = XCAFDoc_DocumentTool_MaterialTool(doc.Main())
    shape_tool = h_shape_tool.GetObject()
    colour_tool = h_colour_tool.GetObject()

    ### Read files
    STEPReader = STEPCAFControl_Reader()
    STEPReader.SetColorMode(True)
#    STEPReader.SetLayerMode(True)
#    STEPReader.SetNameMode(True)
#    STEPReader.SetMatMode(True)
    
    status = STEPReader.ReadFile(filename)
    STEPReader.Reader().PrintCheckLoad(True,0)
    if status == IFSelect_RetDone:
        STEPReader.Transfer(doc.GetHandle())
    
    ### Get root assembly
    shapeLabels = TDF_LabelSequence()
    shape_tool.GetFreeShapes(shapeLabels)
#    print ('Number of shapes at root :%i'%shapeLabels.Length())
    
    ### Here we only have one solid in the model so we directly use the first solid
    ### If there are more, please use the for loop to get all solids
    shape=shape_tool.GetShape(shapeLabels.Value(1))
    
    topo=Topo(shape)
#    for i in range(shapeLabels.Length()):
#        #print(Labels.Value(i+1).Tag())
#        shape=shape_tool.GetShape(shapeLabels.Value(i+1))
#        shapes.append(shape)
    
    return[topo,colour_tool]





def fix_shape(shape):
    """
    """
    sfs=ShapeFix_Shape(shape)
#    sfs.Init(shape)
    sfs.SetPrecision(0.01)
    #sfs.FixEdgeTool
#    sfs.SetMaxTolerance(0.02)
#    w=sfs.FixWireTool().GetObject()
#    l=w.FixDegenerated()
#    t=w.FixSmall(True,0.01)
#    print(t)
##    sfs.SetMinTolerance(0.001)
    sfs.Perform()
    result=sfs.Shape()

#    sfwf=ShapeFix_Wireframe(shape)
#    sfwf.SetModeDropSmallEdges=True
#    print("max tolerance: %f"%sfwf.MaxTolerance())
#    print("min tolerance: %f"%sfwf.MinTolerance())
#    print("Precision: %f"%sfwf.Precision())
#    sfwf.SetPrecision(0.001)
#    sfwf.SetMaxTolerance(0.002)
#    print("max tolerance: %f"%sfwf.MaxTolerance())
#    print("min tolerance: %f"%sfwf.MinTolerance())
#    print("Precision: %f"%sfwf.Precision())
#    sfwf.FixSmallEdges()
#
#
#    result=sfwf.Shape()
    
    return result
def normalize_node_edge_parameter(nodeInfo, topo):
    """
    This replaces the non-normalized parameter of a node on an edge to the normalized parameter.
    
    """
    edges=list(topo.edges())
    print(len(edges))
    for i in range(0,len(edges)):
        (curve,tmin,tmax)=edge_extreme(edges[i])
        
        for lable in list(nodeInfo[2].values())[i]:
            np=(float(nodeInfo[3][lable][3])-tmin)/(tmax-tmin)###calculate normalized parameter
            nodeInfo[3][lable][3]=np


def normalize_node_face_parameter(nodeInfo, topo):
    """
    This replaces the non-normalized parameter of a node on an edge to the normalized parameter.
    
    """
    faces=list(topo.faces())
    for i in range(0,len(faces)):
        [umin,umax,vmin,vmax]=face_extreme(faces[i])
        
        for lable in list(nodeInfo[4].values())[i]:
            normU=(float(nodeInfo[5][lable][3])-umin)/(umax-umin)###calculate normalized parameter
            nodeInfo[5][lable][3]=normU
            normV=(float(nodeInfo[5][lable][4])-vmin)/(vmax-vmin)###calculate normalized parameter
            nodeInfo[5][lable][4]=normV


if __name__ == "__main__":
#    path="E:/3-MyResearchFellowDrive/4-Results/Test/"
    path="E:/3052180/Desktop/"
    
    modelName="Part5"
    shape=read_step_file(path+modelName+".stp")
    topo=Topo(shape)
#    result=fix_shape(shape)
#    wires=topo.wires()
#    for wire in wires:
#        result=fix_shape(wire)
#    print(len(list(topo.edges())))
    fs=list(topo.faces())
    for f in fs:
        vs=list(topo.vertices_from_face(f))
        for v in vs:
            angle=corner_angle(topo,v,f)
            print(angle)
    
    
    