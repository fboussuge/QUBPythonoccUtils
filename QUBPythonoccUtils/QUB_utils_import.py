from OCC.STEPControl import STEPControl_Reader
from OCC.TopAbs import TopAbs_SOLID, TopAbs_FACE, TopAbs_EDGE
from OCC.TopExp import TopExp_Explorer
from OCC.StepRepr import Handle_StepRepr_RepresentationItem

from OCC.TopoDS import topods

import csv

def read_step_file_withnames(filename, breadface = False, breadedge = False):
    """""[Read a step file with names attached to solid, faces (can be extended to edges)]
    
    Arguments:
        filename {[type]} -- [path to the .stp file]
    
    Keyword Arguments:
        breadface {bool} -- [read the faces' names] (default: {False})
        breadedge {bool} -- [read the edges' names] (default: {False})

    Returns:
        (dSolids, dFaces, dEdges) -- [two dicts with name (int) as key and aShape as value]
    """""

    reader = STEPControl_Reader()
    tr = reader.WS().GetObject().TransferReader().GetObject()
    reader.ReadFile(filename)
    reader.TransferRoots()
    shape = reader.OneShape()

    dSolids = dict() #solids initial as a dict with name (int) as key and aShape as value
    dFaces = dict() #faces initial as a dict with name (int) as key and aShape as value
    dEdges = dict() #edges initial as a dict with name (int) as key and aShape as value
    
    #read the solid names
    exp = TopExp_Explorer(shape, TopAbs_SOLID)
    while exp.More():
        s = exp.Current()
        exp.Next()
        #Converting TopoDS_Shape to TopoDS_Face
        solid = topods.Solid(s)
        #Working on the name
        item = tr.EntityFromShapeResult(s, 1)
        if item.IsNull():
            continue
        item = Handle_StepRepr_RepresentationItem.DownCast(item).GetObject()
        name = item.Name().GetObject().ToCString()
        if name:
#            print('Found entity named: {}: {}.'.format(name, s))
            nameid = int(name.split('_')[-1])
            dSolids[nameid] = solid
    
    # read the face names
    if breadface:
        exp = TopExp_Explorer(shape, TopAbs_FACE)
        while exp.More():
            s = exp.Current()
            exp.Next()
            #Converting TopoDS_Shape to TopoDS_Face
            face = topods.Face(s)
            #Working on the name
            item = tr.EntityFromShapeResult(s, 1)
            if item.IsNull():
                continue
            item = Handle_StepRepr_RepresentationItem.DownCast(item).GetObject()
            name = item.Name().GetObject().ToCString()
            if name:
    #            print('Found entity named: {}: {}.'.format(name, s))
                nameid = int(name.split('_')[-1])
                dFaces[nameid] = face

    # read the edge names
    if breadedge:
        exp = TopExp_Explorer(shape, TopAbs_EDGE)
        while exp.More():
            s = exp.Current()
            exp.Next()
            #Converting TopoDS_Shape to TopoDS_Face
            edge = topods.Edge(s)
            #Working on the name
            item = tr.EntityFromShapeResult(s, 1)
            if item.IsNull():
                continue
            item = Handle_StepRepr_RepresentationItem.DownCast(item).GetObject()
            name = item.Name().GetObject().ToCString()
            if name:
    #            print('Found entity named: {}: {}.'.format(name, s))
                nameid = int(name.split('_')[-1])
                dEdges[nameid] = edge

    ret = (dSolids, dFaces, dEdges)
    return ret


def read_csv_file(csvnamefile):
    """[Read a csv file and return the row content]
    
    Arguments:
        csvnamefile {[str]} -- [path of the csv file]
    
    Returns:
        [list] -- [List of int]
    """

     #read the csv file generated from rescal and return the listed entities
    res= []
    with open(csvnamefile, 'r') as csvfile:
       x = csv.reader(csvfile,delimiter=',',quotechar='|')
       for row in x:
           for i in range(0, len(row)-1):
               res.append(int(row[i]))
    return res
