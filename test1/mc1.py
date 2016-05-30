# -*- coding: utf-8 -*-

###
### This file is generated automatically by SALOME v7.7.1 with dump python functionality
###

import sys
import salome

salome.salome_init()
theStudy = salome.myStudy

import salome_notebook
notebook = salome_notebook.NoteBook(theStudy)
sys.path.insert( 0, r'C:/Salome')

###
### GEOM component
###

import GEOM
from salome.geom import geomBuilder
import math
import SALOMEDS


geompy = geomBuilder.New(theStudy)

Vertex_0 = geompy.MakeVertex(0, 0, 0)
Vertex_1 = geompy.MakeVertex(0.2, 0, 0)
Vertex_2 = geompy.MakeVertex(0.8, 0, 0)
Vertex_3 = geompy.MakeVertex(1, 0, 0)
Vertex_4 = geompy.MakeVertex(0, 0.2, 0)
Vertex_5 = geompy.MakeVertex(1, 0.2, 0)
Vertex_6 = geompy.MakeVertex(0, 0.8, 0)
Vertex_7 = geompy.MakeVertex(1, 0.8, 0)
Vertex_8 = geompy.MakeVertex(0, 1, 0)
Vertex_9 = geompy.MakeVertex(0.2, 1, 0)
Vertex_10 = geompy.MakeVertex(0.8, 1, 0)
Vertex_11 = geompy.MakeVertex(1, 1, 0)

Line_0 = geompy.MakeLineTwoPnt(Vertex_0, Vertex_1)
Line_1 = geompy.MakeLineTwoPnt(Vertex_0, Vertex_4)
Line_2 = geompy.MakeLineTwoPnt(Vertex_1, Vertex_2)
Line_3 = geompy.MakeLineTwoPnt(Vertex_1, Vertex_4)

Line_4 = geompy.MakeLineTwoPnt(Vertex_2, Vertex_3)
Line_5 = geompy.MakeLineTwoPnt(Vertex_2, Vertex_5)
Line_6 = geompy.MakeLineTwoPnt(Vertex_3, Vertex_5)
Line_7 = geompy.MakeLineTwoPnt(Vertex_4, Vertex_6)

Line_8 = geompy.MakeLineTwoPnt(Vertex_5, Vertex_7)
Line_9 = geompy.MakeLineTwoPnt(Vertex_6, Vertex_8)
Line_10 = geompy.MakeLineTwoPnt(Vertex_6, Vertex_9)
Line_11 = geompy.MakeLineTwoPnt(Vertex_7, Vertex_10)

Line_12 = geompy.MakeLineTwoPnt(Vertex_7, Vertex_11)
Line_13 = geompy.MakeLineTwoPnt(Vertex_8, Vertex_9)
Line_14 = geompy.MakeLineTwoPnt(Vertex_9, Vertex_10)
Line_15 = geompy.MakeLineTwoPnt(Vertex_10, Vertex_11)

Face1Wire = geompy.MakeWire([Line_0, Line_1, Line_3], 1e-007)
Face2Wire = geompy.MakeWire([Line_2, Line_5, Line_8, Line_11, Line_14, Line_10, Line_7, Line_3], 1e-007)
Face3Wire = geompy.MakeWire([Line_4, Line_5, Line_6], 1e-007)
Face4Wire = geompy.MakeWire([Line_9, Line_10, Line_13], 1e-007)
Face5Wire = geompy.MakeWire([Line_11, Line_12, Line_15], 1e-007)

Face_1 = geompy.MakeFaceWires([Face1Wire], 1)
Face_2 = geompy.MakeFaceWires([Face2Wire], 1)
Face_3 = geompy.MakeFaceWires([Face3Wire], 1)
Face_4 = geompy.MakeFaceWires([Face4Wire], 1)
Face_5 = geompy.MakeFaceWires([Face5Wire], 1)

geompy.addToStudy( Vertex_0, 'Vertex_0' )
geompy.addToStudy( Vertex_1, 'Vertex_1' )
geompy.addToStudy( Vertex_2, 'Vertex_2' )
geompy.addToStudy( Vertex_3, 'Vertex_3' )
geompy.addToStudy( Vertex_4, 'Vertex_4' )
geompy.addToStudy( Vertex_5, 'Vertex_5' )
geompy.addToStudy( Vertex_6, 'Vertex_6' )
geompy.addToStudy( Vertex_7, 'Vertex_7' )
geompy.addToStudy( Vertex_8, 'Vertex_8' )
geompy.addToStudy( Vertex_9, 'Vertex_9' )
geompy.addToStudy( Vertex_10, 'Vertex_10' )
geompy.addToStudy( Vertex_11, 'Vertex_11' )

geompy.addToStudy( Line_0, 'Line_0' )
geompy.addToStudy( Line_1, 'Line_1' )
geompy.addToStudy( Line_2, 'Line_2' )
geompy.addToStudy( Line_3, 'Line_3' )

geompy.addToStudy( Line_4, 'Line_4' )
geompy.addToStudy( Line_5, 'Line_5' )
geompy.addToStudy( Line_6, 'Line_6' )
geompy.addToStudy( Line_7, 'Line_7' )


geompy.addToStudy( Line_8, 'Line_8' )
geompy.addToStudy( Line_9, 'Line_9' )
geompy.addToStudy( Line_10, 'Line_10' )
geompy.addToStudy( Line_11, 'Line_11' )


geompy.addToStudy( Line_12, 'Line_12' )
geompy.addToStudy( Line_13, 'Line_13' )
geompy.addToStudy( Line_14, 'Line_14' )
geompy.addToStudy( Line_15, 'Line_15' )

geompy.addToStudy(Face1Wire, 'Face1Wire')
geompy.addToStudy(Face2Wire, 'Face2Wire')
geompy.addToStudy(Face3Wire, 'Face3Wire')
geompy.addToStudy(Face4Wire, 'Face4Wire')
geompy.addToStudy(Face5Wire, 'Face5Wire')

geompy.addToStudy(Face_1, 'Face_1')
geompy.addToStudy(Face_2, 'Face_2')
geompy.addToStudy(Face_3, 'Face_3')
geompy.addToStudy(Face_4, 'Face_4')
geompy.addToStudy(Face_5, 'Face_5')


###
### SMESH component
###

import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

smesh = smeshBuilder.New(theStudy)

Mesh_1 = smesh.Mesh(Face_1)
NETGEN_2D_1 = Mesh_1.Triangle(algo=smeshBuilder.NETGEN_1D2D)
NETGEN_2D_Parameters_1 = NETGEN_2D_1.Parameters()
NETGEN_2D_Parameters_1.SetMaxSize( 0.04 )
NETGEN_2D_Parameters_1.SetSecondOrder( 0 )
NETGEN_2D_Parameters_1.SetOptimize( 1 )
NETGEN_2D_Parameters_1.SetFineness( 2 )
NETGEN_2D_Parameters_1.SetMinSize( 0 )
NETGEN_2D_Parameters_1.SetUseSurfaceCurvature( 1 )
NETGEN_2D_Parameters_1.SetFuseEdges( 1 )
NETGEN_2D_Parameters_1.SetQuadAllowed( 0 )
isDone = Mesh_1.Compute()


Mesh_2 = smesh.Mesh(Face_2)
NETGEN_2D_2 = Mesh_2.Triangle(algo=smeshBuilder.NETGEN_1D2D)
NETGEN_2D_Parameters_2 = NETGEN_2D_2.Parameters()
NETGEN_2D_Parameters_2.SetMaxSize( 0.04 )
NETGEN_2D_Parameters_2.SetSecondOrder( 0 )
NETGEN_2D_Parameters_2.SetOptimize( 1 )
NETGEN_2D_Parameters_2.SetFineness( 2 )
NETGEN_2D_Parameters_2.SetMinSize( 0 )
NETGEN_2D_Parameters_2.SetUseSurfaceCurvature( 1 )
NETGEN_2D_Parameters_2.SetFuseEdges( 1 )
NETGEN_2D_Parameters_2.SetQuadAllowed( 0 )
Mesh_2.Compute()

Mesh_3 = smesh.Mesh(Face_3)
NETGEN_2D_3 = Mesh_3.Triangle(algo=smeshBuilder.NETGEN_1D2D)
NETGEN_2D_Parameters_3 = NETGEN_2D_3.Parameters()
NETGEN_2D_Parameters_3.SetMaxSize( 0.04 )
NETGEN_2D_Parameters_3.SetSecondOrder( 0 )
NETGEN_2D_Parameters_3.SetOptimize( 1 )
NETGEN_2D_Parameters_3.SetFineness( 2 )
NETGEN_2D_Parameters_3.SetMinSize( 0 )
NETGEN_2D_Parameters_3.SetUseSurfaceCurvature( 1 )
NETGEN_2D_Parameters_3.SetFuseEdges( 1 )
NETGEN_2D_Parameters_3.SetQuadAllowed( 0 )
Mesh_3.Compute()


Mesh_4 = smesh.Mesh(Face_4)
NETGEN_2D_4 = Mesh_4.Triangle(algo=smeshBuilder.NETGEN_1D2D)
NETGEN_2D_Parameters_4 = NETGEN_2D_4.Parameters()
NETGEN_2D_Parameters_4.SetMaxSize( 0.04 )
NETGEN_2D_Parameters_4.SetSecondOrder( 0 )
NETGEN_2D_Parameters_4.SetOptimize( 1 )
NETGEN_2D_Parameters_4.SetFineness( 2 )
NETGEN_2D_Parameters_4.SetMinSize( 0 )
NETGEN_2D_Parameters_4.SetUseSurfaceCurvature( 1 )
NETGEN_2D_Parameters_4.SetFuseEdges( 1 )
NETGEN_2D_Parameters_4.SetQuadAllowed( 0 )
Mesh_4.Compute()

Mesh_5 = smesh.Mesh(Face_5)
NETGEN_2D_5 = Mesh_5.Triangle(algo=smeshBuilder.NETGEN_1D2D)
NETGEN_2D_Parameters_5 = NETGEN_2D_5.Parameters()
NETGEN_2D_Parameters_5.SetMaxSize( 0.04 )
NETGEN_2D_Parameters_5.SetSecondOrder( 0 )
NETGEN_2D_Parameters_5.SetOptimize( 1 )
NETGEN_2D_Parameters_5.SetFineness( 2 )
NETGEN_2D_Parameters_5.SetMinSize( 0 )
NETGEN_2D_Parameters_5.SetUseSurfaceCurvature( 1 )
NETGEN_2D_Parameters_5.SetFuseEdges( 1 )
NETGEN_2D_Parameters_5.SetQuadAllowed( 0 )
Mesh_5.Compute()

FaceGroup_1 = Mesh_1.GroupOnGeom(Face_1,'FaceGroup_1',SMESH.FACE)
FaceGroup_1.SetColor( SALOMEDS.Color( 1, 0.666667, 0 ))
FaceBoundGroup_1 = Mesh_1.GroupOnGeom(Face_1,'FaceBoundGroup_1',SMESH.EDGE)
FaceBoundGroup_1.SetColor( SALOMEDS.Color( 1, 0.666667, 0 ))

FaceGroup_2 = Mesh_2.GroupOnGeom(Face_2,'FaceGroup_2',SMESH.FACE)
FaceGroup_2.SetColor( SALOMEDS.Color( 1, 0.333333, 0 ))
FaceBoundGroup_2 = Mesh_2.GroupOnGeom(Face_2,'FaceBoundGroup_2',SMESH.EDGE)
FaceBoundGroup_2.SetColor( SALOMEDS.Color( 1, 0.333333, 0 ))

FaceGroup_3 = Mesh_3.GroupOnGeom(Face_3,'FaceGroup_3',SMESH.FACE)
FaceGroup_3.SetColor( SALOMEDS.Color( 0.666667, 0.333333, 0 ))
FaceBoundGroup_3 = Mesh_3.GroupOnGeom(Face_3,'FaceBoundGroup_3',SMESH.EDGE)
FaceBoundGroup_3.SetColor( SALOMEDS.Color( 0.666667, 0.333333, 0 ))

FaceGroup_4 = Mesh_4.GroupOnGeom(Face_4,'FaceGroup_4',SMESH.FACE)
FaceGroup_4.SetColor( SALOMEDS.Color( 0.666667, 0.333333, 1 ))
FaceBoundGroup_4 = Mesh_4.GroupOnGeom(Face_4,'FaceBoundGroup_4',SMESH.EDGE)
FaceBoundGroup_4.SetColor( SALOMEDS.Color( 0.666667, 0.333333, 1 ))

FaceGroup_5 = Mesh_5.GroupOnGeom(Face_5,'FaceGroup_5',SMESH.FACE)
FaceGroup_5.SetColor( SALOMEDS.Color( 0.666667, 1, 1 ))
FaceBoundGroup_5 = Mesh_5.GroupOnGeom(Face_5,'FaceBoundGroup_5',SMESH.EDGE)
FaceBoundGroup_5.SetColor( SALOMEDS.Color( 0.666667, 1, 1 ))

try:
	Mesh_1.ExportDAT("FaceGroup_1.dat", FaceGroup_1 )
	Mesh_1.ExportDAT("FaceBoundGroup_1.dat", FaceBoundGroup_1 )
except:
  print 'ExportPartToDAT() of Mesh_1 failed. Invalid file name?'
  
try:
  Mesh_2.ExportDAT("FaceGroup_2.dat", FaceGroup_2 )
  Mesh_2.ExportDAT( "FaceBoundGroup_2.dat", FaceBoundGroup_2 )
except:
  print 'ExportPartToDAT() of Mesh_2 failed. Invalid file name?'
  
try:
  Mesh_3.ExportDAT( "FaceGroup_3.dat", FaceGroup_3 )
  Mesh_3.ExportDAT( "FaceBoundGroup_3.dat", FaceBoundGroup_3 )
except:
  print 'ExportPartToDAT() of Mesh_3 failed. Invalid file name?'
  
try:
  Mesh_4.ExportDAT( "FaceGroup_4.dat", FaceGroup_4 )
  Mesh_4.ExportDAT( "FaceBoundGroup_4.dat", FaceBoundGroup_4 )
except:
  print 'ExportPartToDAT() of Mesh_4 failed. Invalid file name?'
  
try:
  Mesh_5.ExportDAT( "FaceGroup_5.dat", FaceGroup_5 )
  Mesh_5.ExportDAT( "FaceBoundGroup_5.dat", FaceBoundGroup_5 )
except:
  print 'ExportPartToDAT() of Mesh_1 failed. Invalid file name?'

smesh.SetName(Mesh_1.GetMesh(), 'FACE_1_MESH')
smesh.SetName(Mesh_2.GetMesh(), 'FACE_2_MESH')
smesh.SetName(Mesh_3.GetMesh(), 'FACE_3_MESH')
smesh.SetName(Mesh_4.GetMesh(), 'FACE_4_MESH')
smesh.SetName(Mesh_5.GetMesh(), 'FACE_5_MESH')

smesh.SetName(NETGEN_2D_Parameters_1, 'FACE_1_NETGEN_2D_PAR')
smesh.SetName(NETGEN_2D_Parameters_2, 'FACE_2_NETGEN_2D_PAR')
smesh.SetName(NETGEN_2D_Parameters_3, 'FACE_3_NETGEN_2D_PAR')
smesh.SetName(NETGEN_2D_Parameters_4, 'FACE_4_NETGEN_2D_PAR')
smesh.SetName(NETGEN_2D_Parameters_5, 'FACE_5_NETGEN_2D_PAR')

if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser(1)
