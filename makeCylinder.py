#!/usr/bin/env python
# -*- coding: utf-8 -*-

import vtkmodules.all as vtk

def make_cube(x,y,z):

    #define basis vectors
    v1 = [x,0,0]
    v2 = [0,y,0]
    v3 = [0,0,z]

    # Create points for the vertices of the parallelepiped
    points = vtk.vtkPoints()
    points.InsertNextPoint(0.0, 0.0, 0.0)
    points.InsertNextPoint(v1)
    points.InsertNextPoint([v1[i] + v2[i] for i in range(3)])
    points.InsertNextPoint(v2)
    points.InsertNextPoint(v3)
    points.InsertNextPoint([v1[i] + v3[i] for i in range(3)])
    points.InsertNextPoint([v1[i] + v2[i] + v3[i] for i in range(3)])
    points.InsertNextPoint([v2[i] + v3[i] for i in range(3)])

    # Create a VTK cell for the parallelepiped
    parallelepiped = vtk.vtkHexahedron()
    for i in range(8):
        parallelepiped.GetPointIds().SetId(i, i)

    # Create a VTK unstructured grid and add the points and cell to it
    grid = vtk.vtkUnstructuredGrid()
    grid.SetPoints(points)
    grid.InsertNextCell(parallelepiped.GetCellType(), parallelepiped.GetPointIds())

    # Create a mapper
    mapper = vtk.vtkDataSetMapper()
    mapper.SetInputData(grid)

    # Create an actor
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    # Make the faces transparent
    actor.GetProperty().SetOpacity(0.05)

    # Set the color of the edges to black
    actor.GetProperty().SetColor(1.0, 1.0, 1.0)

    # Make the edges thicker
    actor.GetProperty().SetLineWidth(1.0)

    return actor


