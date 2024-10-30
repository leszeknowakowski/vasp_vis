#!/usr/bin/env python
# -*- coding: utf-8 -*-
import time

print('starting...')
tic = time.perf_counter()
import json
toc = time.perf_counter()
print(f'importing json {toc-tic}')

tic = time.perf_counter()
from VASPparser import PoscarParser, OutcarParser
toc = time.perf_counter()
print(f'importing vasp parseer{toc - tic}')

tic = time.perf_counter()
from scipy.spatial.distance import pdist, squareform
toc = time.perf_counter()
print(f'importing scipy {toc - tic}')

tic = time.perf_counter()
import numpy as np
toc = time.perf_counter()
print(f'importing numpy {toc - tic}')

tic = time.perf_counter()
from vtk import *
toc = time.perf_counter()
print(f'importing vtk {toc - tic}')

tic = time.perf_counter()
from makeCylinder import make_cube
toc = time.perf_counter()
print(f'importing cylinder {toc - tic}')

tic = time.perf_counter()
import os
import sys
toc = time.perf_counter()
print(f'importing os & sys {toc - tic}')
tic = time.perf_counter()

# Setting the Qt bindings for QtPy
os.environ["QT_API"] = "PyQt5"
from qtpy import QtWidgets, QtCore, QtGui
from qtpy.QtWidgets import QMainWindow, QApplication, QPushButton, QWidget, QAction, QTabWidget, QVBoxLayout, QLabel, \
    QFileDialog,QLineEdit,QMessageBox, QProgressBar,QDialog
toc = time.perf_counter()
print(f'importing Pyqt {toc - tic}')

tic = time.perf_counter()
from pyvistaqt import QtInteractor, MainWindow
import pyvista as pv
toc = time.perf_counter()
print(f'importing pyvista and pvqt {toc - tic}')

tic = time.perf_counter()
from exceptions import *
toc = time.perf_counter()
print(f'importing exceptions {toc - tic}')

tic = time.perf_counter()
#sys.path.insert(1, "/home/lnowakowski/python/Scripts/")
import chopPARCHG_test_chgcar_comp as chp
toc = time.perf_counter()
print(f'importing CHGCAR parser {toc-tic}')

tic = time.perf_counter()
import matplotlib
matplotlib.use('Qt5Agg')
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure
from matplotlib.colors import ListedColormap
toc = time.perf_counter()
print(f'importing matpplotlib...{toc - tic} s')

from RangeSlider import QRangeSlider
import pathlib


class MplCanvas(FigureCanvasQTAgg):
    def __init__(self, parent=None, width=5, height=7, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = fig.add_subplot(111)
        # self.axes.set_yscale('symlog')
        self.axes.set_xlabel('geometry number')
        self.axes.set_ylabel('energy')
        self.axes.xaxis.tick_top()
        super(MplCanvas, self).__init__(fig)
        
        
class DialogWIndow(QWidget):
    def __init__(self):
        super().__init__()
        layout = QVBoxLayout()
        self.setMinimumSize(300,300)
        self.label1 = QLabel('processing...')
        layout.addWidget(self.label1)
        self.setLayout(layout)   

        
class MyMainWindow(MainWindow):
    """Sets the window and widgets"""

    def __init__(self, parent=None):
                # initialize main window
        tic = time.perf_counter()
        QtWidgets.QMainWindow.__init__(self, parent)
        # create the frame and layout

        self.frame = QtWidgets.QFrame()
        self.frame.setMinimumSize(QtCore.QSize(1280, 720))

        hlayout = QtWidgets.QHBoxLayout(self.frame)
        self.leftLayout = QVBoxLayout()
        hlayout.addLayout(self.leftLayout)

        self.vlayout = QtWidgets.QVBoxLayout()
        self.vlayout.setAlignment(QtCore.Qt.AlignTop)
        self.leftLayout.addLayout(self.vlayout)

        self.vblayout = QtWidgets.QGridLayout()
        self.leftLayout.addLayout(self.vblayout)

        self.tabs = QTabWidget()
        self.tab1 = QWidget()
        self.tab2 = QWidget()
        self.tab3 = QWidget()
        self.tabs.resize(800, 500)

        # Add tabs
        self.tabs.addTab(self.tab2, "Geometry")
        self.tabs.addTab(self.tab1, "")
        self.tabs.addTab(self.tab3, "")

        # Create first tab
        self.tab1.layout = QVBoxLayout(self)


        self.tab1.setLayout(self.tab1.layout)
        self.l = QLabel()
        self.l.setText("This is the first tab")
        self.tab1.layout.addWidget(self.l)
        self.tab1.setLayout(self.tab1.layout)

        # create second tab
        # add the pyvista interactor object
        self.tab2.layout = QVBoxLayout(self)
        self.plotter =QtInteractor()
        self.tab2.layout.addWidget(self.plotter.interactor)
        self.tab2.setLayout(self.tab2.layout)

        # Add tabs to widget
        hlayout.addWidget(self.tabs)

        self.setCentralWidget(self.frame)
        toc = time.perf_counter()
        print(f'initializing main window... {toc - tic}')
        self.show()
    
       
        # initialize variables
        print('initializing application')
        tic = time.perf_counter()
        if not os.path.exists('OUTCAR'):
            print('no OUTCAR found! importing CONTCAR or POSCAR')
            self.outcar_file = False
        else:
            self.outcar_file = True
            self.outcar_data = OutcarParser('OUTCAR')
            self.outcar_coordinates = self.outcar_data.find_coordinates()
            self.outcar_energies = self.outcar_data.find_energy()
        if not os.path.exists('CONTCAR'):
            if not os.path.exists('POSCAR'):
                p = input("eneter file name: ")
                if not os.path.exists(p):
                    raise FileNotFoundError('No important files found! Missing POSCAR')
                else:
                    self.poscar=PoscarParser(p)
                    self.coordinates = self.poscar.coordinates()
                    if not self.outcar_file:
                        self.outcar_coordinates = [self.poscar.coordinates()]
                        self.outcar_energies = [0]

            else:
                self.poscar = PoscarParser('POSCAR')
                self.coordinates = self.poscar.coordinates()
                if not self.outcar_file:
                    self.outcar_coordinates = [self.poscar.coordinates()]
                    self.outcar_energies = [0]
        else:
            if os.path.getsize('CONTCAR') > 0:
                self.poscar = PoscarParser('CONTCAR')
                self.coordinates = self.poscar.coordinates()
                if self.outcar_file == False:
                    self.outcar_coordinates = [self.poscar.coordinates()]
                    self.outcar_energies = [0]

            else:
                if not os.path.exists('POSCAR'):
                    raise EmptyFile('CONTCAR is found but appears to be empty! POSCAR missing! Check your files')
                else:
                    self.poscar = PoscarParser('POSCAR')
                    self.coordinates = self.poscar.coordinates()
                    if not self.outcar_file or self.outcar_coordinates == []:
                        self.outcar_coordinates = [self.poscar.coordinates()]
                        self.outcar_energies = [0]

        self.x = self.poscar.unit_cell_vectors()[0][0]
        self.y = self.poscar.unit_cell_vectors()[1][1]
        self.z = self.poscar.unit_cell_vectors()[2][2]
        self.number_of_atoms = len(self.coordinates)
        self.symbols = self.poscar.list_atomic_symbols()
        self.constrains = self.poscar.constrains()

        script_dir = pathlib.Path(__file__).parent.resolve()
        with open(os.path.join(script_dir, 'elementColorSchemes.json'), 'r') as file:
            color_data = json.load(file)
        self.atom_colors = [color_data[self.symbols[i]] for i in range(self.number_of_atoms)]
        self.bond_threshold = 2.4
        self.coord_pairs = []  # pairs of points connected by a bond
        self.bond_actors = []  # list of bond actors
        self.sphere_actors = []  # list of sphere actors
        self.geometry_actors = []  # list of geometries, each with actors list
        self.constrain_actor = None
        self.planeSource = None
        self.plane_position = int(self.z / 2) * 100
        self.plane_actor = None
        self.plane_actor_heigher = None
        self.symb_actor = None
        self.cube_actor = None
        self.mag_actor = None
        self.bond_actors = None
        self.master_bond_visibility = 2
        self.scatter_item = None
        self.charge_data = None
        self.contour_type = 'total'
        self.eps = 0.1
        self.sphere_radius = 0.5
        toc = time.perf_counter()
        print(f'initializing variables: {toc - tic}')       
        
        ####################################################################### create widgets ####################################################################
        tic = time.perf_counter()
        sphere_cb = QtWidgets.QCheckBox(self.frame)
        sphere_cb.setChecked(True)
        sphere_cb.setText("Sphere")
        
        self.sphere_radius_slider = QtWidgets.QSlider(self.frame)
        self.sphere_radius_slider.setOrientation(QtCore.Qt.Horizontal)
        self.sphere_radius_slider.setMinimum(0)
        self.sphere_radius_slider.setMaximum(10)
        self.sphere_radius_slider.setValue(5)
        self.sphere_radius_slider.setFixedWidth(200)

        numbers_cb = QtWidgets.QCheckBox(self.frame)
        numbers_cb.setChecked(False)
        numbers_cb.setText('numbers')

        self.mag_cb = QtWidgets.QCheckBox(self.frame)
        self.mag_cb.setChecked(False)
        self.mag_cb.setText('magnetization')

        self.constrains_cb = QtWidgets.QCheckBox(self.frame)
        self.constrains_cb.setChecked(False)
        self.constrains_cb.setText('show constrains')

        self.constrains_all_cb = QtWidgets.QCheckBox(self.frame)
        self.constrains_all_cb.setChecked(False)
        self.constrains_all_cb.setText('show all constrains')

        unit_cell_cb = QtWidgets.QCheckBox(self.frame)
        unit_cell_cb.setChecked(True)
        unit_cell_cb.setText('unit cell')

        bonds_cb = QtWidgets.QCheckBox(self.frame)
        bonds_cb.setChecked(True)
        bonds_cb.setText("bonds")

        plane_cb = QtWidgets.QCheckBox(self.frame)
        plane_cb.setChecked(True)
        plane_cb.setText('plane')

        self.plane_height_range_slider = QRangeSlider(self.frame)
        self.plane_height_range_slider.setRange(1, 50)
        self.plane_height_range_slider.setBackgroundStyle('background: qlineargradient(x1:0, y1:0, x2:0, y2:1, stop:0 #222, stop:1 #333);')
        self.plane_height_range_slider.handle.setStyleSheet('background: qlineargradient(x1:0, y1:0, x2:0, y2:1, stop:0 #282, stop:1 #393);')

        self.plane_heigth_slider = QtWidgets.QSlider(self.frame)
        self.plane_heigth_slider.setOrientation(QtCore.Qt.Horizontal)
        self.plane_heigth_slider.setMinimum(0)
        self.plane_heigth_slider.setMaximum(int(self.z) * 100)
        self.plane_heigth_slider.setValue(int(self.z / 2) * 100)
        self.plane_heigth_slider.setFixedWidth(200)

        self.plane_height_label = QtWidgets.QLabel('plane height')

        self.plane_opacity_slider = QtWidgets.QSlider(self.frame)
        self.plane_opacity_slider.setOrientation(QtCore.Qt.Horizontal)
        self.plane_opacity_slider.setMinimum(0)
        self.plane_opacity_slider.setMaximum(100)
        self.plane_opacity_slider.setValue(100)
        self.plane_opacity_slider.setFixedWidth(200)

        self.plane_opacity_label = QtWidgets.QLabel('plane opacity')

        self.geometry_slider = QtWidgets.QSlider(self.frame)
        self.geometry_slider.setOrientation(QtCore.Qt.Horizontal)
        self.geometry_slider.setMinimum(1)
        self.geometry_slider.setMaximum(len(self.outcar_coordinates) - 1)
        self.geometry_slider.setValue(1)
        self.geometry_slider.setTickInterval(1)
        self.geometry_slider.setSingleStep(1)
        self.geometry_slider.setFixedWidth(200)

        self.bond_threshold_slider = QtWidgets.QSlider(self.frame)
        self.bond_threshold_slider.setOrientation(QtCore.Qt.Horizontal)
        self.bond_threshold_slider.setMinimum(200)
        self.bond_threshold_slider.setMaximum(300)
        self.bond_threshold_slider.setValue(int(self.bond_threshold * 100))
        self.bond_threshold_slider.setFixedWidth(200)

        self.bond_threshold_label = QtWidgets.QLabel()
        self.update_bond_threshold_label()

        self.geometry_value_label = QtWidgets.QLabel()
        self.geometry_value_label.setFixedWidth(150)
        self.update_geometry_value_label()
        self.end_geometry_button = QtWidgets.QPushButton()
        self.end_geometry_button.setIcon(QtGui.QIcon('end.png'))
        self.end_geometry_button.setFixedWidth(30)
        self.start_geometry_button = QtWidgets.QPushButton()
        self.start_geometry_button.setIcon(QtGui.QIcon('start.png'))
        self.start_geometry_button.setFixedWidth(30)
        self.next_geometry_button = QtWidgets.QPushButton()
        self.next_geometry_button.setIcon(QtGui.QIcon('next.png'))
        self.next_geometry_button.setFixedWidth(30)
        self.back_geometry_button = QtWidgets.QPushButton()
        self.back_geometry_button.setIcon(QtGui.QIcon('back.png'))
        self.back_geometry_button.setFixedWidth(30)

        self.animate_button = QtWidgets.QPushButton('Animate geometry')
        self.animate_button.setFixedWidth(120)
        self.stop_button = QtWidgets.QPushButton('stop animating')
        self.stop_button.setFixedWidth(120)
        self.stop_button.setEnabled(False)
        
        self.image_seq_button = QtWidgets.QPushButton('image seqence')
        self.image_seq_button.setFixedWidth(200)
        self.make_gif_button = QtWidgets.QPushButton('make gif')
        self.make_gif_button.setFixedWidth(200)
        
        self.energy_plot = MplCanvas(self, width=5, height=4, dpi=100)
        self.energy_plot.axes.plot(self.outcar_energies)

        self.chg_file_button = QPushButton("choose charge density file", self)

        self.chg_eps_slider = QtWidgets.QSlider(self.frame)
        self.chg_eps_slider.setOrientation(QtCore.Qt.Horizontal)
        self.chg_eps_slider.setMinimum(0)
        self.chg_eps_slider.setMaximum(100)
        self.chg_eps_slider.setValue(10)
        self.chg_eps_slider.setTickInterval(1)
        self.chg_eps_slider.setSingleStep(1)
        self.chg_eps_slider.setFixedWidth(800)


        self.total_contours_cb = QPushButton("total density")
        self.spin_contours_cb = QPushButton("spin density")
        self.alfa_contours_cb = QPushButton("alfa density")
        self.beta_contours_cb = QPushButton("beta density")
        self.no_contours_b = QPushButton("clear contours")
        
        ############################################# connect widgets ###############################################################
        sphere_cb.stateChanged.connect(self.toggle_spheres)

        unit_cell_cb.stateChanged.connect(self.toggle_unit_cell)
        self.constrains_cb.stateChanged.connect(self.toggle_constrain_above_plane)
        self.constrains_all_cb.stateChanged.connect(self.toggle_constrain)
        self.plane_heigth_slider.valueChanged.connect(self.toggle_mag_above_plane)
        self.mag_cb.stateChanged.connect(self.toggle_mag_above_plane)
        numbers_cb.stateChanged.connect(self.toggle_symbols)
        bonds_cb.stateChanged.connect(self.toogle_bonds)
        plane_cb.stateChanged.connect(self.toggle_plane)
        self.geometry_slider.valueChanged.connect(self.add_sphere)
        self.geometry_slider.valueChanged.connect(self.update_geometry_value_label)
        self.geometry_slider.valueChanged.connect(self.add_bonds)
        self.geometry_slider.valueChanged.connect(self.add_scatter)
        self.start_geometry_button.clicked.connect(self.start_button)
        self.back_geometry_button.clicked.connect(self.back_geometry)
        self.next_geometry_button.clicked.connect(self.next_geometry)
        self.end_geometry_button.clicked.connect(self.end_geometry)
        self.bond_threshold_slider.valueChanged.connect(self.set_bond_threshold)
        self.bond_threshold_slider.valueChanged.connect(self.add_bonds)
        self.bond_threshold_slider.valueChanged.connect(self.update_bond_threshold_label)
        self.plane_heigth_slider.valueChanged.connect(self.plane_height)
        self.plane_heigth_slider.valueChanged.connect(self.toggle_constrain_above_plane)

        self.plane_height_range_slider.startValueChanged.connect(self.toggle_mag_above_plane)
        self.plane_height_range_slider.startValueChanged.connect(self.all_planes_position)
        self.plane_height_range_slider.startValueChanged.connect(self.toggle_constrain_above_plane)

        self.plane_height_range_slider.endValueChanged.connect(self.toggle_mag_above_plane)
        self.plane_height_range_slider.endValueChanged.connect(self.all_planes_position)
        self.plane_height_range_slider.endValueChanged.connect(self.toggle_constrain_above_plane)

        self.plane_opacity_slider.valueChanged.connect(self.plane_opacity)
        self.animate_button.clicked.connect(self.animate_slider)
        self.stop_button.clicked.connect(self.stop_animation)
        self.image_seq_button.clicked.connect(self.save_image_sequence)
        self.make_gif_button.clicked.connect(self.make_gif)

        self.chg_file_button.clicked.connect(self.show_chg_dialog)

        self.chg_eps_slider.sliderReleased.connect(self.update_eps)
        self.chg_eps_slider.sliderReleased.connect(self.add_contours)

        self.total_contours_cb.clicked.connect(lambda: self.set_contour_type("total"))
        self.total_contours_cb.clicked.connect(self.add_contours)
        self.spin_contours_cb.clicked.connect(lambda: self.set_contour_type("spin"))
        self.spin_contours_cb.clicked.connect(self.add_contours)
        self.alfa_contours_cb.clicked.connect(lambda: self.set_contour_type("alfa"))
        self.alfa_contours_cb.clicked.connect(self.add_contours)
        self.beta_contours_cb.clicked.connect(lambda: self.set_contour_type("beta"))
        self.beta_contours_cb.clicked.connect(self.add_contours)
        self.no_contours_b.clicked.connect(self.clear_contours)
        constrain_layout = QtWidgets.QHBoxLayout()
        constrain_layout.addWidget(self.constrains_cb)
        constrain_layout.addWidget(self.constrains_all_cb)
        constrain_layout.setAlignment(QtCore.Qt.AlignLeft)

        button_layout = QtWidgets.QHBoxLayout()
        button_layout.addWidget(self.animate_button)
        button_layout.addWidget(self.stop_button)
        button_layout.setAlignment(QtCore.Qt.AlignLeft)
        button_layout.setAlignment(QtCore.Qt.AlignLeft)
        
        image_but_layout = QtWidgets.QHBoxLayout()
        image_but_layout.addWidget(self.image_seq_button)
        image_but_layout.addWidget(self.make_gif_button)
        
        self.sphere_radius_slider.valueChanged.connect(self.change_sphere_size)
        self.sphere_radius_slider.valueChanged.connect(self.add_sphere)
        self.sphere_radius_layout = QtWidgets.QHBoxLayout()
        self.sphere_radius_layout.addWidget(self.sphere_radius_slider)
        self.sphere_radius_layout.addWidget(QLabel('sphere radius'))
        self.sphere_radius_layout.setAlignment(QtCore.Qt.AlignLeft)
        
        slider_layout = QtWidgets.QHBoxLayout()
        slider_layout.addWidget(self.geometry_slider)
        slider_layout.addWidget(self.start_geometry_button)
        slider_layout.addWidget(self.back_geometry_button)
        slider_layout.addWidget(self.next_geometry_button)
        slider_layout.addWidget(self.end_geometry_button)
        slider_layout.addWidget(self.geometry_value_label)
        slider_layout.setAlignment(QtCore.Qt.AlignLeft)

        bond_slider_layout = QtWidgets.QHBoxLayout()
        bond_slider_layout.addWidget(self.bond_threshold_slider)
        bond_slider_layout.addWidget(self.bond_threshold_label)

        plane_height_layout = QtWidgets.QHBoxLayout()
        plane_height_layout.addWidget(self.plane_heigth_slider)
        plane_height_layout.addWidget(self.plane_height_label)

        plane_opacity_layout = QtWidgets.QHBoxLayout()
        plane_opacity_layout.addWidget(self.plane_opacity_slider)
        plane_opacity_layout.addWidget(self.plane_opacity_label)

        chg_epsilon_layout = QtWidgets.QHBoxLayout()
        chg_epsilon_layout.addWidget(QLabel("epsilon: "))
        chg_epsilon_layout.addWidget(self.chg_eps_slider)

        ############################################################################ add widgets ##########################################
        self.vlayout.addWidget(sphere_cb)
        self.vlayout.addLayout(self.sphere_radius_layout)
        self.vlayout.addWidget(numbers_cb)
        self.vlayout.addWidget(unit_cell_cb)
        self.vlayout.addWidget(self.mag_cb)
        self.vlayout.addLayout(constrain_layout)
        self.vlayout.addWidget(bonds_cb)
        self.vlayout.addWidget(plane_cb)
        self.vlayout.addLayout(plane_height_layout)
        self.vlayout.addWidget(self.plane_height_range_slider)
        self.vlayout.addLayout(plane_opacity_layout)
        self.vlayout.addLayout(bond_slider_layout)
        self.vlayout.addLayout(slider_layout)
        self.vlayout.addLayout(button_layout)
        self.vlayout.addLayout(image_but_layout)
        self.vlayout.addWidget(self.add_scatter())

        self.vblayout.addWidget(self.chg_file_button,0,0)
        self.vblayout.addLayout(chg_epsilon_layout,1,0)
        self.vblayout.addWidget(self.total_contours_cb, 2,0)
        self.vblayout.addWidget(self.spin_contours_cb, 2,1)
        self.vblayout.addWidget(self.alfa_contours_cb, 3,0)
        self.vblayout.addWidget(self.beta_contours_cb, 3,1)
        self.vblayout.addWidget(self.no_contours_b, 4,0)
        
        self.timer = QtCore.QTimer()
        self.timer.setInterval(50)
        self.timer.timeout.connect(self.increment_slider)
        toc = time.perf_counter()
        print(f'initializing widgets... {toc - tic}')

        tic = time.perf_counter()
        self.add_unit_cell()
        toc = time.perf_counter()
        print(f'adding unit cell.. {toc - tic}')

        tic = time.perf_counter()
        self.add_sphere()
        toc = time.perf_counter()
        print(f'adding spheres... {toc - tic}')

        tic = time.perf_counter()
        self.add_constrains()
        toc = time.perf_counter()
        print(f'adding constrains... {toc - tic}')

        tic = time.perf_counter()
        self.add_symbol_and_number()
        toc = time.perf_counter()
        print(f'adding symbols... {toc - tic}')

        tic = time.perf_counter()
        self.add_bonds()
        toc = time.perf_counter()
        print(f'adding bonds...{toc - tic} s')

        self.add_plane(self.z / 2)
        self.add_plane_higher(self.z)

        self.plotter.add_key_event('z', self.camera_z)
        self.plotter.add_key_event('x', self.camera_x)
        self.plotter.add_key_event('y', self.camera_y)

        tic = time.perf_counter()
        
        toc = time.perf_counter()
        print(f'show application... {toc - tic}')

        self.plotter.view_yz()
        
    def create_chgcar_data(self):
        tic = time.perf_counter()       
        chopping_factor = 1
        self.charge_data = chp.PoscarParser(self.chg_file_path, chopping_factor)
        self.add_contours()
        
        toc = time.perf_counter()
        print(f'charge density data read and displayed. Time elapsed: {toc - tic} s')        
        
    def set_contour_type(self, contour_type):
        self.contour_type = contour_type

    def update_eps(self):
        self.eps = self.chg_eps_slider.value()/100

    def add_contours(self):
        if self.charge_data == None:
            print("no data was found")
        else:
            if self.contour_type == "total":
                test_vol = self.charge_data.all_numbers[0]
            elif self.contour_type == "spin":
                test_vol = self.charge_data.all_numbers[1]
            elif self.contour_type == "alfa":
                test_vol = self.charge_data.alfa
            elif self.contour_type == "beta":
                test_vol = self.charge_data.beta
            test_vol = test_vol.swapaxes(0, 2)
            max_val = np.max(test_vol)
            min_val = np.min(test_vol)
            pvgrid = pv.ImageData()
            pvgrid.dimensions = test_vol.shape
            pvgrid.origin = (0, 0, 0)
            pvgrid.spacing = tuple(self.charge_data.calculate_grid_spacings())
            pvgrid.point_data["values"] = test_vol.flatten(order="F")
            self.end_geometry()
            self.plane_actor.GetProperty().SetOpacity(0)
            if self.contour_type == "spin":
                contours = pvgrid.contour([self.eps * max_val, self.eps * min_val])
            else:
                contours = pvgrid.contour([self.eps*max_val])
                
            colors = pv.LookupTable()
            colors.scalar_range = (self.eps * min_val, self.eps * max_val)
            lightblue = np.array([0/256,255/256,254/256, 1.0])
            yellow = np.array([255/256,255/256,0/256, 1.0])
            mapping = np.linspace(self.eps * min_val, self.eps * max_val, 256)
            newcolors = np.empty((256,4))
            newcolors[mapping > 0] = yellow
            newcolors[mapping < 0 ] = lightblue
            my_colormap =  ListedColormap(newcolors)
            #colors.cmap = [(0,255,254), 'yellow']
            self.plotter.add_mesh(contours, name='isosurface', smooth_shading=True, opacity=1, cmap=my_colormap)

    def clear_contours(self):
        pvgrid = pv.ImageData()
        pvgrid.dimensions=(10,10,10)
        pvgrid.origin=(0,0,0)
        pvgrid.spacing=(1,1,1)
        data = np.random.rand(10,10,10)
        pvgrid.point_data["values"]=data.flatten(order="F")
        contours = pvgrid.contour()
        self.plotter.add_mesh(contours, opacity=0, name='isosurface')
        
    def show_chg_dialog(self):
        self.w = DialogWIndow()
        self.w.setWindowFlags(QtCore.Qt.WindowStaysOnTopHint)
        self.w.show()
        
        file_dialog = QFileDialog()
        file_path, _ = file_dialog.getOpenFileName(self, "choose charge density file")
        self.chg_file_path = file_path
        
        
        self.create_chgcar_data()
        self.w.close()
    
    def save_image(self):
        self.plotter.screenshot(str(self.geometry_slider.value())+'.png')
        
    def save_image_sequence(self):
        """loops the geometry slider when ends"""
        value = 0
        while value < self.geometry_slider.maximum():
            self.geometry_slider.setValue(value)
            self.save_image()
            value += 1
    def make_gif(self):
        import glob
        from PIL import Image
        png_files = sorted(glob.glob('*png'), key=lambda x: int(x.split('.')[0]))
        frames = []
        for i in png_files:
            new=Image.open(i)
            frames.append(new)
        frames[0].save('geometry.gif', format='GIF', append_images = frames[1:], save_all=True, duration=200, loop=0)
        for i in png_files:
            os.remove(i)       
        
        
    def camera_z(self):
        self.plotter.view_xy()

    def camera_x(self):
        self.plotter.view_yz()

    def camera_y(self):
        self.plotter.view_xz()

    def add_scatter(self):
        """
        removes oryginal graph widget and plots a new one, with same line plot
        and a scatter with acutal energy value referenced to visible geometry.
        Takes ~0.02 s
        """
        self.vlayout.removeWidget(self.energy_plot)
        val = self.geometry_slider.value()
        self.energy_plot = MplCanvas(self, width=3, height=3, dpi=100)
        self.energy_plot.axes.plot(range(1, len(self.outcar_energies)+1),self.outcar_energies)
        self.energy_plot.axes.scatter(val, self.outcar_energies[val-1])
        self.vlayout.addWidget(self.energy_plot)
        return self.energy_plot

    def animate_slider(self):
        """Sets the timer of animation and enables or disabled animation buttons"""
        self.animate_button.setEnabled(False)
        self.stop_button.setEnabled(True)
        self.timer.start()

    def increment_slider(self):
        """loops the geometry slider when ends"""
        value = self.geometry_slider.value()
        value += 1
        if value > self.geometry_slider.maximum():
            value = self.geometry_slider.minimum()
        self.geometry_slider.setValue(value)
        

    def stop_animation(self):
        """stops the animation by stopping the timer"""
        self.timer.stop()
        self.animate_button.setEnabled(True)
        self.stop_button.setEnabled(False)

    def update_geometry_value_label(self):
        """updates the label indicating current geometry number"""
        self.geometry_value_label.setText(f'geometry number: {self.geometry_slider.value()}')

    def update_bond_threshold_label(self):
        """ updates the label indicating current bond visibility value"""
        self.bond_threshold_label.setText(f'bond visibility: {self.bond_threshold_slider.value() / 100}')

    def start_button(self):
        self.geometry_slider.setValue(0)

    def back_geometry(self):
        value = self.geometry_slider.value()
        value -= 1
        self.geometry_slider.setValue(value)

    def next_geometry(self):
        value = self.geometry_slider.value()
        value += 1
        self.geometry_slider.setValue(value)

    def end_geometry(self):
        last = len(self.outcar_coordinates)
        self.geometry_slider.setValue(last)
    def change_sphere_size(self, value):
            self.sphere_radius = value/10
            
    def add_sphere(self):
        """adds atoms from single geometry as spheres to renderer.
         Using VTK code because it is 100x faster than pyvista
         """
        for actor in self.sphere_actors:
            self.plotter.renderer.RemoveActor(actor)
        self.sphere_actors = []  # make local later
        coordinates = self.outcar_coordinates[self.geometry_slider.value()]

        for coord, col in zip(coordinates, self.atom_colors):
            sphere = vtkSphereSource()
            sphere.SetRadius(self.sphere_radius)
            sphere.SetThetaResolution(30)
            sphere.SetPhiResolution(30)
            sphere.SetCenter(*coord)

            sphere_mapper = vtkPolyDataMapper()
            sphere_mapper.SetInputConnection(sphere.GetOutputPort())
            sphere_actor = vtkActor()
            sphere_actor.SetMapper(sphere_mapper)
            sphere_actor.GetProperty().SetColor(col[0] / 255, col[1] / 255, col[2] / 255)
            self.sphere_actors.append(sphere_actor)
            self.plotter.renderer.AddActor(sphere_actor)
        return self.sphere_actors

    def slide_geometry(self, value):
        """
        toggle spheres on and off depending on what geometry is to be seen. First, switches off all
        spheres visibility, and then swtiches on only these, which are pointed by
        geometry_slide.value() (here called as value)
        """
        for geometry in self.geometry_actors:
            for actor in geometry:
                actor.SetVisibility(False)
        for actor in self.geometry_actors[value]:
            actor.SetVisibility(True)

    def add_bonds(self):
        """
        render bonds as cylinder. First calculate all pairs, which distance is less than value,
        and then uses module make_cylinder to render a cylinder between points
        """

        counter = []
        rendercounter = []
        coordinates = self.outcar_coordinates[self.geometry_slider.value()]
        if self.bond_actors is not None:
            for actor in self.bond_actors:
                self.plotter.renderer.RemoveActor(actor)

        self.bond_actors = []
        self.coord_pairs = []
        bond_threshold = self.bond_threshold

        # Calculate pairwise distances
        distances = squareform(pdist(coordinates))

        # Get the upper triangular part of the distances matrix
        upper_triangle = np.triu(distances, k=1)

        # Find pairs with distance less than threshold (excluding distances between the same point pairs)
        pairs = np.argwhere((upper_triangle < bond_threshold) & (upper_triangle > 0))
        colors = vtkNamedColors()
        for pair in pairs:
            point1, point2 = pair
            coord1 = coordinates[point1]
            coord2 = coordinates[point2]
            self.coord_pairs.append([coord1, coord2])
        for pair in self.coord_pairs:
            if self.master_bond_visibility == 2:
                line_source = vtkLineSource()
                line_source.SetPoint1(pair[0])
                line_source.SetPoint2(pair[1])
                mapper = vtkPolyDataMapper()
                mapper.SetInputConnection(line_source.GetOutputPort())
                actor = vtkActor()
                actor.SetMapper(mapper)
                actor.GetProperty().SetLineWidth(5)
                actor.GetProperty().SetColor(colors.GetColor3d('Black'))

                self.plotter.renderer.AddActor(actor)
                self.bond_actors.append(actor)

    def set_bond_threshold(self, value):
        """ setter of the bond_threshold_value"""
        self.bond_threshold = value / 100

    def add_plane(self, value):
        """renders a plane perpendicular to XY plane at value height"""
        if self.plane_actor is not None:
            self.plotter.remove_actor(self.plane_actor)
        colors = vtkNamedColors()
        colors.SetColor('BkgColor', [26, 51, 77, 255])
        self.planeSource = vtkPlaneSource()
        self.planeSource.SetNormal(0.0, 0.0, 1.0)
        self.planeSource.SetOrigin(-5, -5, value)
        self.planeSource.SetPoint1(self.x + 5, -5, value)
        self.planeSource.SetPoint2(-5, self.y + 5, value)
        self.planeSource.Update()
        plane = self.planeSource.GetOutput()

        # Create a mapper and actor
        mapper = vtkPolyDataMapper()
        mapper.SetInputData(plane)
        self.plane_actor = vtkActor()
        self.plane_actor.SetMapper(mapper)
        self.plane_actor.GetProperty().SetColor(colors.GetColor3d('Banana'))
        #  self.plane_actor.GetProperty().SetOpacity()

        self.plotter.renderer.AddActor(self.plane_actor)

    def add_plane_higher(self, value):
        """renders a plane perpendicular to XY plane at value height"""
        if self.plane_actor_heigher is not None:
            self.plotter.remove_actor(self.plane_actor_heigher)
        colors = vtkNamedColors()
        colors.SetColor('BkgColor', [26, 51, 77, 255])
        self.planeSource_heigher = vtkPlaneSource()
        self.planeSource_heigher.SetNormal(0.0, 0.0, 1.0)
        self.planeSource_heigher.SetOrigin(-5, -5, value)
        self.planeSource_heigher.SetPoint1(self.x + 5, -5, value)
        self.planeSource_heigher.SetPoint2(-5, self.y + 5, value)
        self.planeSource_heigher.Update()
        plane = self.planeSource_heigher.GetOutput()

        # Create a mapper and actor
        mapper = vtkPolyDataMapper()
        mapper.SetInputData(plane)
        self.plane_actor_heigher = vtkActor()
        self.plane_actor_heigher.SetMapper(mapper)
        self.plane_actor_heigher.GetProperty().SetColor(colors.GetColor3d('Banana'))
        #  self.plane_actor_heigher.GetProperty().SetOpacity()

        self.plotter.renderer.AddActor(self.plane_actor_heigher)

    def add_unit_cell(self):
        """renders an parallelpipe representig an unit cell"""
        self.cube_actor = make_cube(self.x, self.y, self.z)
        self.cube_actor.visibility = True
        self.plotter.renderer.AddActor(self.cube_actor)
        #self.chgplotter.renderer.AddActor(self.cube_actor)

    def add_constrains(self):
        """ renders a constrains of atoms"""
        self.plotter.renderer.RemoveActor(self.constrain_actor)
        coordinates = self.outcar_coordinates[self.geometry_slider.value()]
        self.constrain_actor = self.plotter.add_point_labels(coordinates, self.constrains, font_size=30,
                                                             show_points=False, always_visible=True, shape=None)
        self.constrain_actor.SetVisibility(False)

    def add_symbol_and_number(self):
        """ renders an atom symbol and number"""
        self.plotter.renderer.RemoveActor(self.symb_actor)
        coordinates = self.outcar_coordinates[self.geometry_slider.value()]
        symb_num = self.poscar.symbol_and_number()
        coords = list(np.array(coordinates) + [0, -0.2, 0.5])
        self.symb_actor = self.plotter.add_point_labels(coords, symb_num, font_size=30, show_points=False,
                                                        always_visible=False, shape=None)
        self.symb_actor.SetVisibility(False)

    def toggle_spheres(self, flag):
        """switches on and off spheres visibility"""
        for actor in self.sphere_actors:
            actor.SetVisibility(flag)

    def toggle_unit_cell(self, flag):
        """ switches on and off unit cell visibility"""
        self.cube_actor.SetVisibility(flag)

    def toggle_constrain(self, flag):
        """ switches on and off constrains visibility"""
        self.plotter.renderer.RemoveActor(self.constrain_actor)
        self.add_constrains()
        self.constrain_actor.SetVisibility(flag)

    def toggle_constrain_above_plane(self, flag):
        self.plotter.renderer.RemoveActor(self.constrain_actor)
        if self.constrains_cb.isChecked():
            slidervalue = self.plane_height_range_slider.getRange()
            height = slidervalue[0] / 100 * self.z
            end = slidervalue[1] / 100 * self.z
            coordinates = np.array(self.outcar_coordinates[self.geometry_slider.value()])
            coords = []
            constr = []
            indice = np.where((coordinates[:, 2] > height) & (coordinates[:,2] < end))[0]
            indices = indice.tolist()
            for i in range(len(indices)):
                coords.append(list(coordinates[indices[i]]))
                constr.append(self.constrains[indices[i]])
            self.constrain_actor = self.plotter.add_point_labels(coords, constr, font_size=30,
                                                                 show_points=False, always_visible=True, shape=None)
            self.constrain_actor.SetVisibility(flag)

    def all_planes_position(self, value):
        val = self.plane_height_range_slider.getRange()
        startVal = val[0]
        endVal = val[1]
        self.plane_position = startVal
        self.planeSource.SetOrigin(-5, -5, startVal / 100 * self.z)
        self.planeSource.SetPoint1(self.x + 5, -5, startVal / 100 * self.z)
        self.planeSource.SetPoint2(-5, self.y + 5, startVal / 100 * self.z)
        self.planeSource.Update()

        self.planeSource_heigher.SetOrigin(-5, -5, endVal / 100 * self.z)
        self.planeSource_heigher.SetPoint1(self.x + 5, -5, endVal / 100 * self.z)
        self.planeSource_heigher.SetPoint2(-5, self.y + 5, endVal / 100 * self.z)
        self.planeSource_heigher.Update()

        self.plotter.renderer.Render()

    def toggle_mag_above_plane(self, flag):
        self.plotter.renderer.RemoveActor(self.mag_actor)
        if self.mag_cb.isChecked():
            self.end_geometry()
            slidervalue = self.plane_height_range_slider.getRange()
            height = slidervalue[0] / 100 * self.z
            end = slidervalue[1] / 100 * self.z

            coordinates = np.array(self.outcar_coordinates[self.geometry_slider.value()])
            coords = []
            magnet = []
            mag = self.outcar_data.find_magnetization()
            indices = list(np.where((coordinates[:, 2] > height) & (coordinates[:,2] < end))[0])
            for i in range(len(indices)):
                coords.append(list(coordinates[indices[i]]))
                magnet.append(mag[indices[i]])
            self.mag_actor = self.plotter.add_point_labels(coords, magnet, font_size=30,
                                                                 show_points=False, always_visible=True, shape=None)
            self.mag_actor.SetVisibility(flag)

    def toggle_symbols(self, flag):
        """ switches on and off symbols and numbers visibility"""
        self.symb_actor.SetVisibility(flag)

    def toogle_bonds(self, flag):
        """ switches on and off bonds visibility"""
        for actor in self.bond_actors:
            actor.SetVisibility(flag)
        self.master_bond_visibility = flag

    def plane_height(self, value):
        """ changes plane height """
        self.plane_position = value
        self.planeSource.SetOrigin(-5, -5, value / 100)
        self.planeSource.SetPoint1(self.x + 5, -5, value / 100)
        self.planeSource.SetPoint2(-5, self.y + 5, value / 100)
        self.planeSource.Update()
        self.plotter.renderer.Render()

    def plane_opacity(self, value):
        """ changes plane opacity"""
        self.plane_actor.GetProperty().SetOpacity(value / 100)

    def toggle_plane(self, flag):
        """ switches on and off plane visibility"""
        self.plane_actor.GetProperty().SetOpacity(flag)


if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)
    window = MyMainWindow()
    sys.exit(app.exec_())
