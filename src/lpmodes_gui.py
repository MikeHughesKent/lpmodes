# -*- coding: utf-8 -*-
"""
lpmodes_gui

A Graphical User Interface (GUI) for interacting with the lpmodes package. The
GUI is based on PyQT5. Modes are displayed using the ImageDisplayQT widget.

Parameters of the fibre (NAs, core radius) and wavelength are defined in the 
top left. Clicking 'Find Modes' then finds solutions. The table in the centre
is updated with the mode details. Power in core is not calculated immediately 
as this is more time consuming, click 'Find Power in Core' button to update
this. Select a line of the table to display an image of that mode on the 
right. There are options on the left to choose the image grid size, the 
physical size of the image as a multiple of the core size, whether to display
the amplitude or the intensity, and whether to display a circle showing the 
core size. The data can be exported as CSV, or the current mode image or all
mode images saved as tifs using buttons on the left.

@author: Michael Hughes, Applied Optics Group, University of Kent
"""

import sys 
import os
import inspect
from pathlib import Path
import time
import math

import numpy as np
import matplotlib.pyplot as plt

from PyQt5 import QtGui, QtCore, QtWidgets   
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtWidgets import *
from PyQt5.QtGui import QIcon, QPalette, QColor, QImage, QPixmap, QPainter, QPen, QGuiApplication
from PyQt5.QtGui import QPainter, QBrush, QPen
from PyQt5.QtXml import QDomDocument, QDomElement

from PIL import Image, TiffImagePlugin

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), 'gui_lib')))

from image_display import ImageDisplay
from im_tools import  save_image8, save_tif_stack

import lpmodes

class LPModesGui(QMainWindow):
    
    def __init__(self, parent = None):
        
        super(LPModesGui, self).__init__(parent)
        
        self.create_layout()
        self.setMinimumWidth(1200)
        self.clear()        
        
        
    def create_layout(self):
        """ Set up QT
        """
        
        self.setWindowTitle("LP Modes")
        
        self.main_widget = QWidget()
        self.setCentralWidget(self.main_widget)
        self.main_widget.setLayout(main_layout:=QHBoxLayout())
        
        self.menu_panel = QWidget()
        self.menu_panel.setMinimumWidth(300)
        self.menu_panel.setMaximumWidth(300)

        vert = QVBoxLayout()
        self.menu_panel.setLayout(vert)
        
        ##### Physical Parameters    
        title = QLabel("Mode display Options:")
        font = title.font()
        font.setBold(True)
        title.setFont(font)
        vert.addWidget(title)
        grid = QGridLayout()

        vert.addLayout(grid)
                        
        self.wavelength_spin = QDoubleSpinBox()
        self.wavelength_spin.setMaximum(5000)
        
        self.core_n_spin = QDoubleSpinBox()
        self.core_n_spin.setMinimum(1)
        self.core_n_spin.setMaximum(2)
        
        self.cladding_n_spin = QDoubleSpinBox()
        self.cladding_n_spin.setMinimum(1)
        self.cladding_n_spin.setMaximum(2)
        
        self.core_radius_spin = QDoubleSpinBox()
        self.core_radius_spin.setMaximum(1000)
        
        self.wavelength_spin.setValue(0.5)
        self.core_n_spin.setValue(1.4)
        self.cladding_n_spin.setValue(1.38)
        self.core_radius_spin.setValue(5)
                
        grid.addWidget(QLabel("Wavelength (um):"),0,0)
        grid.addWidget(QLabel("Core n:"),1,0)
        grid.addWidget(QLabel("Cladding n:"), 2,0)
        grid.addWidget(QLabel("Core Radius (um):"),3,0)
        
        grid.addWidget(self.wavelength_spin,0,1)
        grid.addWidget(self.core_n_spin,1,1)
        grid.addWidget(self.cladding_n_spin,2,1)
        grid.addWidget(self.core_radius_spin,3,1)
        
        self.find_modes_btn = QPushButton("Find Modes")        
        self.find_modes_btn.clicked.connect(self.find_modes_clicked)
        
        self.find_power_btn = QPushButton("Find Power in Modes")        
        self.find_power_btn.clicked.connect(self.find_power_clicked)
     
        #### Report Modes 
        grid2 = QGridLayout()
        grid2.addWidget(QLabel("Numerical Aperture:"), 0, 0)
        grid2.addWidget(QLabel("Predicted Modes:"), 1, 0)
        grid2.addWidget(QLabel("Modes Found:"), 2, 0)
        grid2.addWidget(QLabel("Rot. Modes Found:"), 3, 0)
        
        self.na_disp = QLineEdit()
        self.na_disp.setEnabled(False)
        #self.na_disp.setMaximumWidth(100)
        grid2.addWidget(self.na_disp, 0, 1)    

        
        self.predicted_modes_disp = QLineEdit()
        self.predicted_modes_disp.setEnabled(False)
        grid2.addWidget(self.predicted_modes_disp, 1, 1)
        
        self.found_modes_disp = QLineEdit()
        self.found_modes_disp.setEnabled(False)
        grid2.addWidget(self.found_modes_disp, 2, 1)
        
        self.rot_modes_disp = QLineEdit()
        self.rot_modes_disp.setEnabled(False)
        grid2.addWidget(self.rot_modes_disp, 3, 1)
             
        self.info_box = QLabel("")
        self.info_box.setWordWrap(True)
        self.info_box.setStyleSheet("border: 1px solid black;")

        vert.addWidget(self.find_modes_btn)
        vert.addWidget(self.find_power_btn)
        vert.addItem(QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding))

        vert.addLayout(grid2)

        vert.addItem(QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding))

        #### Mode Display Options
        title2 = QLabel("Mode display Options:")
        font = title2.font()
        font.setBold(True)
        title2.setFont(font)
        vert.addWidget(title2)

        self.display_amplitude_check = QCheckBox("Display Amplitude")
        self.display_amplitude_check.stateChanged.connect(self.change_mode_display)
        self.show_core_check = QCheckBox("Show Core Diameter")
        self.show_core_check.stateChanged.connect(self.change_mode_display)
        self.show_rotated_check = QCheckBox("Show Rotation")
        self.show_rotated_check.stateChanged.connect(self.change_mode_display)

                                         
        grid3 = QGridLayout()
        grid3.addWidget(QLabel("Grid Size (px):"),0,0)
        
        self.grid_size_spin = QSpinBox()
        self.grid_size_spin.setMaximum(1000)
        self.grid_size_spin.setMinimum(10)
        self.grid_size_spin.setValue(100)
        self.grid_size_spin.valueChanged[int].connect(self.change_mode_display)


        self.size_factor_spin = QDoubleSpinBox()
        self.size_factor_spin.setMaximum(10)
        self.size_factor_spin.setMinimum(1)
        self.size_factor_spin.setValue(1.5)
        self.size_factor_spin.valueChanged[float].connect(self.change_mode_display)
        
        grid3.addWidget(QLabel("Grid Size (px):"),0,0)
        grid3.addWidget(self.grid_size_spin, 0, 1)
        grid3.addWidget(QLabel("Display Size (x radius):"),1,0)
        grid3.addWidget(self.size_factor_spin, 1, 1)
        
        vert.addLayout(grid3)
        vert.addWidget(self.display_amplitude_check)
        vert.addWidget(self.show_core_check)
        vert.addWidget(self.show_rotated_check)


        #### Export Buttons
        vert.addItem(QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding))

        self.export_data_btn = QPushButton("Export data as CSV")
        self.export_data_btn.clicked.connect(self.export_data)
        vert.addWidget(self.export_data_btn)
        
        self.export_mode_btn = QPushButton("Save current mode image")
        self.export_mode_btn.clicked.connect(self.export_mode)
        vert.addWidget(self.export_mode_btn)

        self.export_stack_btn = QPushButton("Save stack of mode images")
        self.export_stack_btn.clicked.connect(self.export_stack)
        vert.addWidget(self.export_stack_btn)

        
        #### Table
        self.tableWidget = QTableWidget() 
        self.tableWidget.setMinimumWidth(400)
        self.tableWidget.setRowCount(0)  
        self.tableWidget.setColumnCount(6) 
        
        self.tableWidget.setSelectionBehavior(QTableView.SelectRows)
        selection_model = self.tableWidget.selectionModel()
        selection_model.selectionChanged.connect(self.table_selection_changed)
        self.tableWidget.setSelectionMode(QTableView.SingleSelection)
          
        #Table will fit the screen horizontally 
        self.tableWidget.horizontalHeader().setStretchLastSection(True) 
        self.tableWidget.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch) 
        self.tableWidget.setStyleSheet("""
                                       QTableWidget {
                                           padding-right: 10px;  # Add 10 pixels of space on the right side
                                           }
                                       """)
        self.tableWidget.setHorizontalHeaderLabels([
            "l", "m", "u", "Beta", "n eff", "Core P"
        ])
        
        self.tableWidget.setColumnWidth(0, 20)  # Width of the second column ("l")
        self.tableWidget.setColumnWidth(1, 20)  # Width of the third column ("m")
        self.tableWidget.setColumnWidth(2, 30)  # Width of the fourth column ("u")
        self.tableWidget.setColumnWidth(3, 30) # Width of the fifth column ("Beta")
        self.tableWidget.setColumnWidth(4, 30) # Width of the sixth column ("n eff")
        self.tableWidget.setColumnWidth(5, 50) # Width of the sixth column ("n eff")


        table_layout = QHBoxLayout()

        table_layout.addWidget(self.tableWidget)

        author = QLabel("LPModes GUI 1.0, Michael Hughes, Applied Optics Group, Physics & Astronomy, University of Kent.")
        author.setWordWrap(True)
        vert.addItem(QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding))

        vert.addStretch()
        vert.addWidget(author)        
 
        #### Mode Display
        displayFrame = QWidget()
        displayFrame.setLayout(layout:=QVBoxLayout())
        self.image_display = ImageDisplay()
        self.image_display.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)

        layout.addWidget(self.image_display)
        policy = QSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        displayFrame.setSizePolicy(policy)


        #### Final assembly
        main_layout.addWidget(self.menu_panel, stretch = 0)
        
        # Create a QWidget to hold the table and spacer
        table_widget_container = QWidget()
        table_widget_container.setLayout(table_layout)        
             
        splitter = QSplitter()
        splitter.addWidget(table_widget_container)
        splitter.addWidget(displayFrame)

        splitter.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        splitter.setStretchFactor(0, 1)  # Table has a lower stretch factor (less priority to shrink)
        splitter.setStretchFactor(1, 2)  # Image display has a higher stretch factor (more priority to expand)
        splitter.setCollapsible(0, False)
        splitter.setCollapsible(1, False)
        splitter.setHandleWidth(2)  # Increase the width of the splitter handle        
        splitter.setStyleSheet("""
        QSplitter::handle {
        background-color: #A0A0A0; /* Darker gray handle */
        }
        """)

        main_layout.addWidget(self.menu_panel)  # The menu panel remains outside the splitter
        main_layout.addWidget(splitter)  # Add the splitter to the main layout
       
 
    def find_modes_clicked(self):
        """ Finds all allowed modes and updates table
        """
        
        self.clear()
        
        self.core_radius = self.core_radius_spin.value()
        self.core_n = self.core_n_spin.value()
        self.cladding_n = self.cladding_n_spin.value()
        self.wavelength = self.wavelength_spin.value()
        
        QApplication.setOverrideCursor(Qt.WaitCursor)
        self.modes = lpmodes.find_modes(self.core_radius, self.core_n, self.cladding_n, self.wavelength)
        QApplication.restoreOverrideCursor()

        #### Mode info updates
        num_modes = len(self.modes)
        num_modes_rot = lpmodes.num_rotated_modes(self.modes)
        est_modes = round(lpmodes.est_num_modes(self.core_radius, self.core_n, self.cladding_n, self.wavelength) / 2)
        na = round(lpmodes.fibre_na(self.core_n, self.cladding_n),2)
        self.na_disp.setText(f"{na}")
        self.predicted_modes_disp.setText(f"{est_modes}")
        self.found_modes_disp.setText(f"{num_modes}")
        self.rot_modes_disp.setText(f"{num_modes_rot}")

        self.change_mode_display()

        #### Table update
        self.tableWidget.setRowCount(len(self.modes))  

        for idx,mode in enumerate(self.modes):
            self.tableWidget.setItem(idx,0, QTableWidgetItem(str(mode.l))) 
            self.tableWidget.setItem(idx,1, QTableWidgetItem(str(mode.m))) 
            self.tableWidget.setItem(idx,2, QTableWidgetItem(str(round(mode.u,3))))
            self.tableWidget.setItem(idx,3, QTableWidgetItem(str(round(mode.beta,3))))
            self.tableWidget.setItem(idx,4, QTableWidgetItem(str(round(mode.n_eff,3))))
        self.tableWidget.selectRow(0)
        self.tableWidget.setFocus()


    def find_power_clicked(self):
        """ Calculates power in core for each mode and adds to table
        """
        
        QApplication.setOverrideCursor(Qt.WaitCursor)
        
        for idx,mode in enumerate(self.modes):
            power_in_core = lpmodes.power_in_core(mode, self.core_radius * self.size_factor_spin.value(), self.grid_size_spin.value())
            self.tableWidget.setItem(idx, 5, QTableWidgetItem(str(round(power_in_core,3))))
        QApplication.restoreOverrideCursor()
        self.tableWidget.selectRow(0)
        self.tableWidget.setFocus()
            
            
    def change_mode_display(self):
        """ Draws the current selected mode
        """
        
        self.image_display.clear_overlays()

        if self.modes is not None:
        
            #if self.solution is not None:
            if self.show_rotated_check.isChecked():
                 im = self.modes[self.mode_selected].plot_amplitude_rotated(self.core_radius * self.size_factor_spin.value(), self.grid_size_spin.value())
            else:
                 im = self.modes[self.mode_selected].plot_amplitude(self.core_radius * self.size_factor_spin.value(), self.grid_size_spin.value())

            if self.display_amplitude_check.isChecked():
                 self.image_display.set_colormap(lpmodes.ampcol())
                 self.image_display.autoScale = False
                 im = im / np.max(np.abs(im)) * 127
                 im = im + 127
            else:
                 self.image_display.set_colormap('hot')
                 self.image_display.autoScale = False
                 im = np.abs(im)**2
                 im = im/ np.max(im) * 255

            self.image_display.set_image(im)

            if self.show_core_check.isChecked():
                pixel_size = self.grid_size_spin.value() / (self.core_radius * self.size_factor_spin.value())
                circle_size = self.core_radius * pixel_size
                loc = (self.grid_size_spin.value() - circle_size) / 2
                overlay = self.image_display.add_overlay(ImageDisplay.ELLIPSE, loc, loc, circle_size, circle_size, QPen(Qt.blue, 2, Qt.SolidLine), QBrush())
       
      
    def table_selection_changed(self, line):
        """ Handles a change in which line of the table is selected
        """
        
        mode_selected = line.indexes()[0].row() 
        if mode_selected >= 0 and mode_selected < len(self.modes):
            self.mode_selected = mode_selected
        self.change_mode_display()        
        
                
    def clear(self):
        """ Wipes table and mode image
        """
      
        self.info_box.setText("Modes not calculated")
        self.solution = None
        self.modes = None
        self.mode_selected = 0
        self.image_display.set_image(np.zeros((100,100)))
        
            
    def export_data(self):
        """ Writes mode data to CSV
        """
        
        if self.modes is not None:
            try:
                filename = QFileDialog.getSaveFileName(self, 'Select filename to save to:', '', filter='*.csv')[0]
            except:
                filename = None 
            if filename is not None:    
                lpmodes.modes_to_csv(self.modes, filename)
            else:
                QMessageBox.about(self, "Error", "Cannot write to file.") 
                

    def export_mode(self):
        """ Saves image of currently selected mode
        """
        
        if self.modes is not None:
            try:
                filename = QFileDialog.getSaveFileName(self, 'Select filename to save to:', '', filter='*.tif')[0]
            except:
                filename = None 
                
            if filename is not None:  
                if self.show_rotated_check.isChecked():
                    im = self.modes[self.mode_selected].plot_amplitude_rotated(self.core_radius * self.size_factor_spin.value(), self.grid_size_spin.value())
                else:
                    im = self.modes[self.mode_selected].plot_amplitude(self.core_radius * self.size_factor_spin.value(), self.grid_size_spin.value())
                    
                if self.display_amplitude_check.isChecked():
                    im = im / np.max(np.abs(im)) * 127
                    im = im + 127 
                else:
                    im = np.abs(im)**2
                    im = im / np.max(im) * 255     
                save_image8(im, filename)
            else:
               QMessageBox.about(self, "Error", "Cannot save image.") 
            
            
    
    def export_stack(self):
        """ Saves images of all modes to a tif stack
        """
        
        if self.modes is not None:
            try:
                filename = QFileDialog.getSaveFileName(self, 'Select filename to save to:', '', filter='*.tif')[0]
            except:
                filename = None 
                
            if filename is not None: 
                if self.show_rotated_check.isChecked():                    
                    stack = lpmodes.plot_modes_amplitude(self.modes, self.core_radius * self.size_factor_spin.value(), self.grid_size_spin.value())[1]
                else:
                    stack = lpmodes.plot_modes_amplitude(self.modes, self.core_radius * self.size_factor_spin.value(), self.grid_size_spin.value())[0]
                    
                if self.display_amplitude_check.isChecked():
                    stack = stack / np.max(np.abs(stack)) * 127
                    stack = stack + 127 
                else:
                    stack = np.abs(stack)**2
                    stack = stack / np.max(stack) * 255 
                    
                QApplication.setOverrideCursor(Qt.WaitCursor)
                time.sleep(0.2)
                save_tif_stack(stack, filename)
                QApplication.restoreOverrideCursor()

            else:
               QMessageBox.about(self, "Error", "Cannot save images.") 
               

    
# Launch the GUI
if __name__ == '__main__':
        
   app=QApplication(sys.argv)
   app.setStyle("Fusion")
     
   # Create and display GUI
   window = LPModesGui()
   window.show()
   
   # When the window is closed, close everything
   sys.exit(app.exec_())

    
