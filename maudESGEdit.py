#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Make it work in python 2 or 3
from __future__ import print_function

"""
Copyright (C) S. Merkel, Universite de Lille, France

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
"""

# Some text for the GUI. I put it on top so it is easier to upate
bottomWindowLabel = "Utility to fix data in esg files for MAUD (c) 2020-now, S. Merkel, Univ. Lille"
aboutWindowText = """
<h3>MAUD ESG editor</h3>
(c) 2020-now, S. Merkel, Univ. Lille
<P>Utility to fix data in ESG files before Rietveld refinement in MAUD.</P>
<P>This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.</P>
<P>The source code is available online at https://github.com/smerkel/maudESGEdit/.</P>
"""

helpWindowText = """
<h3>MAUD ESG editor</h3>
(c) 2020, S. Merkel, Univ. Lille
<P>This programs allows for<ul>
<li>Removing rubbish data points, due to detector gaps, for instance,</li>
<li>Removing and correction your background.</li>
</ul>

<h4>Rubbish data points</h4>
<P>Rubbish data points occur when you  have either a parasite signal or gaps in your data collection. MAUD can exclude 2theta ranges, but will fail if the parasite signals occur at inconsistent 2theta values. This utility allows to browse through all azimuth in an ESG file and remove these rubbish data points.</P>
<P>How to proceed to remove rubbish data points?<ul>
<li>Click on the zoom icon before you start, to activate the zoom mode,</li>
<li>Zoom-in on the data points you wish to remove,</li>
<li>Click on <i>Remove data points</i> or hit <i>Ctr-d</i>.</li>
</ul>

<h4>Background</h4>
<P>First, if your background can be treated in MAUD, treat them in MAUD. It is always a million times better to fit your signal within a single software with clear hypothesis for the modeling of experimental data. This utility is for hopeless cases, inconsistent detector chips inducing jumps in the overall background level, for instance.</P>
<P>How to subtract background?<ul>
<li>Right click to generate background points (background interpolation is linear),</li>
<li>Select <i>Subtract background</i> when you have enough points.</li>
</ul>

<h4>Baseline intensity</h4>
<P>MAUD will fail with negative or zero intensity points. It actually is a feature, but can really throw you off with some datas. Your refinement will diverge. To fix this, You can<ul>
<li>Add a fixed intensity to the data at the current azimuth by selecting <i>Shift dataset...</i>,</li>
<li>Add a fixed intensity to the data at all azimuth by selecting <i>Shift all datasets...</i>,</li>
<li>Set the minimum intensity at all azimuth by selecting <i>Set value for minimum intensity...</i>. For each azimuth, the program will evaluate the current intensity minimum and add the necessary shift to bring the minimum intensity to the value decided by the user. This operation is typically performed at the end, once you have edited data for all azimuth.</li>
</ul>

<h4>Navigation and file management</h4>
<P>You can navigate between spectra using the <i>Previous</i> or <i>Next</i> buttons or by using the <i>left</i> and <i>right</i> keyboard keys.</P>
<P>When you are done, save your new data in a new ESG file.</P>

<h4>Masks</h4>
<P>As you remove rubbish data points, we record 2theta ranges along with the corresponding azimuths. You can save these ranges in a file, with a <i>msk</i> extension to reuse them later.</P>
<P>If you want to remove the same data ranges as in a previous processing, use the <i>Mask -> Load and apply mask<i>menu item.</P>

<h4>Final note</h4>
<P>Is this data manipulation? If you use this sotware to remove actual data, it is. If you use this software to clean up spectra (due to gaps in your detectors, for instance), it is not.</P>
<P>Good luck with your data!</P>
"""

# System functions, to manipulate command line arguments
import sys
import argparse
from argparse import RawTextHelpFormatter
import os.path

# Plotting routines
import matplotlib
matplotlib.use("Qt5Agg")
from matplotlib.figure import Figure
from matplotlib.backend_bases import key_press_handler, Event
from matplotlib.backends.backend_qt5agg import (
    FigureCanvasQTAgg as FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar)

# PyQT graphical interface
import PyQt5.QtWidgets 
import PyQt5.QtCore
import PyQt5.QtGui
import base64

# Maths stuff
import numpy
import math
import scipy
from scipy import interpolate

# Baseline removal tools. Removed it. Does not work with our drops in intensity
#from skued import baseline_dt

# Useful stuff
import copy

#def RunningMedian(x,N):
    #idx = numpy.arange(N) + numpy.arange(len(x)-N+1)[:,None]
    #b = [row[row>0] for row in x[idx]]
    #return numpy.array(map(numpy.median,b))

#################################################################
#
# Save a read MAUD esg files, mask files
#
#################################################################


def parseESG(filename):
	# Reads number of spectra found
	text = open(filename).read()
	nspectra = int(text.count('loop_'))
	# print ("Looking at file %s I found %d spectra." % (filename, nspectra))
	# Read the ESG file, Reads all the lines and saves it to the array "content"
	f = open(filename, 'r')
	logcontent = [line.strip() for line in f.readlines()]
	f.close()
	# Locating lines with new spectra. They start with _pd_block_id. We look for lines with this
	lookup = "_pd_block_id"
	linesdb = []
	for num, line in enumerate(logcontent, 0):
		if lookup in line:
			linesdb.append(num)
	# For each spectrum, save header, etaangle, and data
	headers = []
	data = []
	etas = []
	detdistance = 0.
	i = 0
	for linestart in linesdb:
		header = ""
		line = linestart
		test = True
		while test:
			txt = logcontent[line]
			header = header + txt + "\n"
			if (txt == "_pd_meas_intensity_total"):
				test = False # We reached the end headers
			else:
				# We search for information we need
				a=txt.split()
				if (len(a) > 0):
					if (a[0] == "_pd_instr_dist_spec/detc"):
						detdistance = float(a[1])
					if (a[0] == "_pd_meas_angle_eta"):
						eta = a[1]
			line += 1
		headers.append(header)
		etas.append(eta)
		# Done with header for this spectrum, now look at data
		test = True
		thisdata = []
		while test:
			txt = logcontent[line]
			a=txt.split()
			if (len(a) < 2):
				test = False
			else:
				twotetha = math.degrees(math.atan(float(a[0])/detdistance))
				thisdata.append([twotetha, float(a[0]), float(a[1])])
			line += 1
		#print ("Read data %d, found %d lines with intensities" % (i, len(thisdata)))
		i += 1
		data.append(thisdata)
	toreturn = {}
	toreturn["headers"] = headers
	toreturn["data"] = data
	toreturn["etas"] = etas
	toreturn["detdistance"] = detdistance
	return toreturn

def saveEsgToFile(esgData,filename):
	headers = esgData["headers"]
	data = esgData["data"]
	neta = len(headers)
	string = ""
	for i in range(0,neta):
		string += headers[i]
		thisdata = data[i]
		for j in range(0,len(thisdata)):
			if (not(numpy.isnan(thisdata[j][2]))):
				string += "%.2f %.8f\n" % (thisdata[j][1], thisdata[j][2])
		string += "\n"
	# Ready to save
	f = open(filename, 'w')
	f.write(string)
	f.close()
	return

def saveMaskToFile(mask, filename):
	string = "# Version: 1.0\n# Mask for maudESGEdit\n# Each line: azimuth number, 2theta range to remove\n# Looks a bit like a cif file, but not a true CIF\n#\n"
	string += "\nloop_\n_esg_azimuth_number _esg_2theta_delete_min _esg_2theta_delete_max\n"
	items = []
	# Creation command was self.mask.append({"set":True, "eta": self.etaToPlot, "clear2thetamin": left, "clear2thetamax": right})
	# Removing empty items
	for item in mask:
		if (item["set"]):
			items.append(item)
	# Sorting mask
	items = sorted(items, key=lambda d: (d['eta'],d['clear2thetamin'])) 
	# Adding to string
	for item in items:
		string += "%i %f %f\n" % (item["eta"], item["clear2thetamin"], item["clear2thetamax"])
	string += "\n"
	# Ready to save
	f = open(filename, 'w')
	f.write(string)
	f.close()
	return

def loadMaskFromFile(filename):
	mask = []
	# Read the Mask file, Reads all the lines and saves it to the array "content"
	f = open(filename, 'r')
	logcontent = [line.strip() for line in f.readlines()]
	f.close()
	# Locating mask data. They start with _esg_azimuth_number. We look for lines with this
	lookup = "_esg_azimuth_number"
	linesdb = []
	for num, line in enumerate(logcontent, 0):
		if lookup in line:
			linesdb.append(num)
	for linestart in linesdb:
		line = linestart + 1 # Starting 3 lines below marker
		test = True
		while test:
			elts = (logcontent[line]).split()
			if (len(elts) < 3):
				test = False # We reached the end
			else:
				# We search for information we need
				mask.append({"set":True, "eta": int(elts[0]), "clear2thetamin": float(elts[1]), "clear2thetamax": float(elts[2])})
			line += 1
	return mask
	

#################################################################
#
# Simple text window
#
#################################################################

class textWindow(PyQt5.QtWidgets.QDialog):
	def __init__(self, title, content, parent=None):
		super(textWindow, self).__init__(parent)
		
		self.main_frame = PyQt5.QtWidgets.QWidget()
		self.setWindowTitle(title) 

		# Add text field
		self.b = PyQt5.QtWidgets.QTextEdit(self)
		self.b.setHtml(content)
		
		# connect button to function on_click
		button = PyQt5.QtWidgets.QPushButton("Close window", self)
		button.clicked.connect(self.on_click)
		hlay = PyQt5.QtWidgets.QHBoxLayout()
		hlay.addItem(PyQt5.QtWidgets.QSpacerItem(300, 10, PyQt5.QtWidgets.QSizePolicy.Expanding))
		hlay.addWidget(button)
		
		# Vertical box layout in the window and setting it up
		vbox = PyQt5.QtWidgets.QVBoxLayout()
		vbox.addWidget(self.b)
		#self.button.setAlignment(PyQt5.QtCore.Qt.AlignRight | PyQt5.QtCore.Qt.AlignVCenter)
		vbox.addLayout(hlay)
		
		# We are done...
		self.setLayout(vbox)
		self.setGeometry(300, 200, 600, 400)
		self.show()

	def on_click(self):
		self.close()

#################################################################
#
# Special dialog to input a twotetha range
#
#################################################################

class TThethaRangeDialog(PyQt5.QtWidgets.QDialog):
	def __init__(self, parent=None):
		super(TThethaRangeDialog, self).__init__(parent)
		self.setWindowTitle("Restrict 2theta range")
		
		self.ok = False

		self.start = PyQt5.QtWidgets.QLineEdit(self)
		# self.start.setValidator(PyQt5.QtGui.QDoubleValidator()) getting lost between French and English, let's stick to English numbers
		self.end = PyQt5.QtWidgets.QLineEdit(self)
		buttonBox = PyQt5.QtWidgets.QDialogButtonBox(PyQt5.QtWidgets.QDialogButtonBox.Ok | PyQt5.QtWidgets.QDialogButtonBox.Cancel, self);

		layout = PyQt5.QtWidgets.QFormLayout(self)
		layout.addRow("Minimum value for 2theta", self.start)
		layout.addRow("Maximum value for 2theta", self.end)
		layout.addWidget(buttonBox)

		buttonBox.accepted.connect(self.accept)
		buttonBox.rejected.connect(self.reject)
		self.show()

	def accept(self):
		try:
			start = float(self.start.text())
		except Exception:
			PyQt5.QtWidgets.QMessageBox.critical(self, 'Error','Minimum value for 2theta is not a number')
			return
		try:
			end = float(self.end.text())
		except Exception:
			PyQt5.QtWidgets.QMessageBox.critical(self, 'Error','Maximum value for 2theta is not a number')
			return
		if (end < start):
			PyQt5.QtWidgets.QMessageBox.critical(self, 'Error','Not a proper range')
			return
		self.ok = True
		self.close()

	def getInputs(self):
		return (float(self.start.text()), float(self.end.text()))
	
	def isOk(self):
		return self.ok


#################################################################
#
# Class to build the Graphical User Interface
#
#################################################################

# Various matplotlib tricks to adapt the GUI to what we want
#
# Access the forward and backward keys in mathplotlib and use them to move between grains
# Inspired from 
# - https://stackoverflow.com/questions/14896580/matplotlib-hooking-in-to-home-back-forward-button-events
# - https://stackoverflow.com/questions/37506260/adding-an-item-in-matplotlib%C2%B4s-toolbar
#
# Click on peak and get information on h,k,l, diffraction angles, and indexing errors
#

#def new_forward(self, *args, **kwargs):
	#s = 'forward_event'
	#event = Event(s, self)
	#event.foo = 100
	#self.canvas.callbacks.process(s, event)
	## forward(self, *args, **kwargs) # If you wanted to still call the old forward event

#def new_backward(self, *args, **kwargs):
	#s = 'backward_event'
	#event = Event(s, self)
	#event.foo = 100
	#self.canvas.callbacks.process(s, event)
	## backward(self, *args, **kwargs) # If you wanted to still call the old backward event

def new_home(self, *args, **kwargs):
    s = 'home_event'
    event = Event(s, self)
    event.foo = 100
    self.canvas.callbacks.process(s, event)
    #home(self, *args, **kwargs)

NavigationToolbar.toolitems = (
	('Home', 'Reset original view', 'home', 'home'), 
	#('Back', 'Previous spectrum', 'back', 'back'), 
	#('Forward', 'Next spectrum', 'forward', 'forward'), 
	#(None, None, None, None), 
	('Pan', 'Pan axes with left mouse, zoom with right', 'move', 'pan'), 
	('Zoom', 'Zoom to rectangle', 'zoom_to_rect', 'zoom'), 
	#(None, None, None, None), 
	('Save', 'Save the figure', 'filesave', 'save_figure'))
	#(None, None, None, None))

NavigationToolbar.home = new_home

class plotEsg(PyQt5.QtWidgets.QMainWindow):
	
	"""
	Constructor
	
	Parmeters:
	
	"""
	def __init__(self, parent=None):
		# Prepare the main window
		PyQt5.QtWidgets.QMainWindow.__init__(self, parent)
		pm = PyQt5.QtGui.QPixmap()
		pm.loadFromData(base64.b64decode(iconXPMbase64))
		i = PyQt5.QtGui.QIcon()
		i.addPixmap(pm)
		self.setWindowIcon(PyQt5.QtGui.QIcon(i))
		# Setting starting variables
		self.nEta = 0				# Number of azimuths
		self.title = "MAUD ESG edit" # Window title
		self.etaToPlot = 0			# Which azimuth are we looking at
		self.needToSave = False		# Set True when something has been changed in the data
		self.olddata = []			# For cancel actions, so we can go back
		self.oldetaToPlot = []		# For cancel actions, so we can go back
		self.mask = []				# Saving mask, 2 theta ranges and azimuth at which to remove data
		self.fileSaveHint = None	# Hint for file saving
		self.xbg = []				# Used for creating background
		self.ybg = []				# Used for creating background
		self.dounzoom = True		# Unzoom when replotting default)
		self.nautobg = 5			# Number of points for auto-background
		self.doautobg = False		# Shall we do autobg?
		self.pathtomask = None		# Path no mask file
		# Done setting variables, preparing the gui
		self.create_main_frame()
		self.on_draw()
		self.show()
	"""
	Builds up the GUI
	"""
	def create_main_frame(self):
		
		# Preparing a frame
		self.main_frame = PyQt5.QtWidgets.QWidget()
		self.setWindowTitle(self.title)
		
		# Building a menu bar
		
		mainMenu = self.menuBar()
		fileMenu = mainMenu.addMenu('File')
		editMenu = mainMenu.addMenu('Edit data')
		bgMenu = mainMenu.addMenu('Background')
		maskMenu = mainMenu.addMenu('Mask')
		helpMenu = mainMenu.addMenu('Help')
		
		openButton = PyQt5.QtWidgets.QAction(PyQt5.QtGui.QIcon.fromTheme("document-open"), 'Open ESG...', self)
		openButton.setShortcut('Ctrl+O')
		openButton.setStatusTip('Open MAUD ESG file...')
		openButton.triggered.connect(self.open_esg)
		fileMenu.addAction(openButton)
		
		saveButton = PyQt5.QtWidgets.QAction(PyQt5.QtGui.QIcon.fromTheme("document-save-as"), 'Save new ESG...', self)
		saveButton.setShortcut('Ctrl+S')
		saveButton.setStatusTip('Save new data for processing in MAUD')
		saveButton.triggered.connect(self.save_esg)
		fileMenu.addAction(saveButton)
		
		fileMenu.addSeparator()
		
		exitButton = PyQt5.QtWidgets.QAction(PyQt5.QtGui.QIcon.fromTheme("application-exit"), 'Exit', self)
		exitButton.setShortcut('Ctrl+Q')
		exitButton.setStatusTip('I am done!')
		exitButton.triggered.connect(self.closeEvent)
		fileMenu.addAction(exitButton)
		
		self.delPointsButton = PyQt5.QtWidgets.QAction(PyQt5.QtGui.QIcon.fromTheme("edit-delete"), 'Remove data points', self)
		self.delPointsButton.setShortcut('Ctrl+D')
		self.delPointsButton.setStatusTip('This will remove all data points within the plot below. Zoom into the region where you want the points removed.')
		self.delPointsButton.triggered.connect(self.remove_points)
		self.delPointsButton.setDisabled(True)
		editMenu.addAction(self.delPointsButton)
		
		self.tthetaButton = PyQt5.QtWidgets.QAction(PyQt5.QtGui.QIcon.fromTheme("zoom-original"), 'Restrict 2theta range', self)
		self.tthetaButton.setShortcut('Ctrl+T')
		self.tthetaButton.setStatusTip('Restrict data to a give 2theta range.')
		self.tthetaButton.triggered.connect(self.edit_twothetarange)
		self.tthetaButton.setDisabled(True)
		editMenu.addAction(self.tthetaButton)
		
		self.cancelButton = PyQt5.QtWidgets.QAction(PyQt5.QtGui.QIcon.fromTheme("edit-undo"), 'Undo', self)
		self.cancelButton.setShortcut('Ctrl+Z')
		self.cancelButton.setStatusTip('I screwed up!')
		self.cancelButton.triggered.connect(self.cancel_last)
		self.cancelButton.setDisabled(True)
		editMenu.addAction(self.cancelButton)
		
		#redoButton = PyQt5.QtWidgets.QAction(PyQt5.QtGui.QIcon.fromTheme("edit-redo"), 'Redo', self)
		#redoButton.setShortcut('Ctrl+R')
		#redoButton.setStatusTip('I did not screw up!')
		#redoButton.triggered.connect(self.redo_last)
		#editMenu.addAction(redoButton)
		
		delButton = PyQt5.QtWidgets.QAction(PyQt5.QtGui.QIcon.fromTheme("edit-delete"), 'Remove background points', self)
		delButton.setStatusTip('This will remove the background points.')
		delButton.triggered.connect(self.remove_bg_points)
		bgMenu.addAction(delButton)
		
		self.subtractBgButton = PyQt5.QtWidgets.QAction(PyQt5.QtGui.QIcon.fromTheme("edit-cut"), 'Subtract background', self)
		self.subtractBgButton.setStatusTip('This will subtract the linear background in blue from the data. Right-click on the mouse to add background points.')
		self.subtractBgButton.setShortcut('Ctrl+B')
		self.subtractBgButton.triggered.connect(self.subtract_background)
		self.subtractBgButton.setDisabled(True)
		bgMenu.addAction(self.subtractBgButton)
		
		delButton = PyQt5.QtWidgets.QAction(PyQt5.QtGui.QIcon.fromTheme("go-up"), 'Shift dataset...', self)
		delButton.setStatusTip('This will add a given value to all intensities for the current dataset.')
		delButton.triggered.connect(self.shift_data)
		bgMenu.addAction(delButton)
		
		delButton = PyQt5.QtWidgets.QAction(PyQt5.QtGui.QIcon.fromTheme("go-top"), 'Shift all datasets...', self)
		delButton.setStatusTip('This will add a given value to all intensities for all datasets.')
		delButton.triggered.connect(self.shift_data_all)
		bgMenu.addAction(delButton)
		
		delButton = PyQt5.QtWidgets.QAction(PyQt5.QtGui.QIcon.fromTheme("go-bottom"), 'Set minimum intensity...', self)
		delButton.setStatusTip('This will set the minimum intensity value for all datasets.')
		delButton.triggered.connect(self.setmin_data_all)
		bgMenu.addAction(delButton)
		
		loadMaskButton = PyQt5.QtWidgets.QAction(PyQt5.QtGui.QIcon.fromTheme("document-open"), 'Load and apply mask...', self)
		loadMaskButton.setShortcut('Ctrl+M')
		loadMaskButton.setStatusTip('Load and apply a mask')
		loadMaskButton.triggered.connect(self.load_and_apply_mask)
		maskMenu.addAction(loadMaskButton)
		
		saveMaskButton = PyQt5.QtWidgets.QAction(PyQt5.QtGui.QIcon.fromTheme("document-save-as"), 'Save mask...', self)
		saveMaskButton.setStatusTip('Save the mask to reuse later')
		saveMaskButton.triggered.connect(self.save_mask)
		maskMenu.addAction(saveMaskButton)
		
		aboutButton = PyQt5.QtWidgets.QAction(PyQt5.QtGui.QIcon.fromTheme("help-contents"), 'User manual...', self)
		aboutButton.setShortcut('Ctrl+H')
		aboutButton.setStatusTip('What is this thing?!')
		aboutButton.triggered.connect(self.helpwindow)
		helpMenu.addAction(aboutButton)
		
		aboutButton = PyQt5.QtWidgets.QAction(PyQt5.QtGui.QIcon.fromTheme("help-about"), 'About this program...', self)
		aboutButton.setShortcut('Ctrl+A')
		aboutButton.setStatusTip('What is this thing?!')
		aboutButton.triggered.connect(self.about)
		helpMenu.addAction(aboutButton)
		
		# Creating a matplotlibe figure

		self.fig = Figure((8.0, 8.0), dpi=100,tight_layout=True,edgecolor='w',facecolor='w')
		
		# Top of the gui with buttons and labels

		self.etaLabel = PyQt5.QtWidgets.QLabel("Spectrum (0-%d) : " % (self.nEta-1),self)
		self.etaNBox = PyQt5.QtWidgets.QLineEdit("%d" % (self.etaToPlot), self)
		self.etaNBox.returnPressed.connect(self.new_eta)
		buttonP = PyQt5.QtWidgets.QPushButton('Previous', self)
		buttonP.setToolTip('Move to previous spectrum')
		buttonP.clicked.connect(self.handle_backward)
		buttonN = PyQt5.QtWidgets.QPushButton('Next', self)
		buttonN.setToolTip('Move to next spectrum')
		buttonN.clicked.connect(self.handle_forward)
		# deleteLabel = PyQt5.QtWidgets.QLabel("Remove data points", self)
		buttonD = PyQt5.QtWidgets.QPushButton("Remove data points", self)
		buttonD.setToolTip('This will remove all data points within the plot below. Zoom into the region where you want the points removed.')
		buttonD.clicked.connect(self.remove_points)
		buttonBg = PyQt5.QtWidgets.QPushButton("Subtract background", self)
		buttonBg.setToolTip('This will subtract the linear background in blue from the data. Right-click on the mouse to add background points.')
		buttonBg.clicked.connect(self.subtract_background)
		hlay = PyQt5.QtWidgets.QHBoxLayout()
		hlay.addWidget(self.etaLabel)
		hlay.addWidget(buttonP)
		hlay.addWidget(self.etaNBox)
		hlay.addWidget(buttonN)
		hlay.addStretch(1)
		#hlay.addItem(PyQt5.QtWidgets.QSpacerItem(300, 10, PyQt5.QtWidgets.QSizePolicy.Expanding))
		#hlay.addWidget(deleteLabel)
		hlay.addWidget(buttonD)
		hlay.addWidget(buttonBg)
		
		# Horizontal layout for autobackground option
		hlay2 = PyQt5.QtWidgets.QHBoxLayout()
		self.autobgbox = PyQt5.QtWidgets.QCheckBox("Auto-background",self)
		self.autobgbox.stateChanged.connect(self.changeautobg)
		lab = PyQt5.QtWidgets.QLabel("Number of auto-bg points",self)
		self.autobgNBox = PyQt5.QtWidgets.QLineEdit("%d" % (self.nautobg), self)
		self.autobgNBox.returnPressed.connect(self.changenautobg)
		if (not self.doautobg):
			self.autobgNBox.setDisabled(True)
		hlay2.addWidget(self.autobgbox)
		hlay2.addWidget(lab)
		hlay2.addWidget(self.autobgNBox)
		hlay2.addStretch(1)
		
		
		
		# Adding keyboard shortcuts
		#shortcut1 = PyQt5.QtWidgets.QShortcut(PyQt5.QtGui.QKeySequence("Ctrl+d"), self)
		#shortcut1.activated.connect(self.remove_points)
		shortcut2 = PyQt5.QtWidgets.QShortcut(PyQt5.QtGui.QKeySequence("Ctrl+n"), self)
		shortcut2.activated.connect(self.handle_forward)
		shortcut3 = PyQt5.QtWidgets.QShortcut(PyQt5.QtGui.QKeySequence("Ctrl+p"), self)
		shortcut3.activated.connect(self.handle_backward)
		shortcut4 = PyQt5.QtWidgets.QShortcut(PyQt5.QtCore.Qt.Key_Right, self)
		shortcut4.activated.connect(self.handle_forward)
		shortcut5 = PyQt5.QtWidgets.QShortcut(PyQt5.QtCore.Qt.Key_Left, self)
		shortcut5.activated.connect(self.handle_backward)
		
		# Adding a canvas for the plot
		self.canvas = FigureCanvas(self.fig)
		self.canvas.setParent(self.main_frame)
		self.canvas.setFocusPolicy(PyQt5.QtCore.Qt.StrongFocus)
		self.canvas.mpl_connect('button_press_event', self.on_press) 
		self.canvas.setFocus()

		# Adding a toolbar and trying to deal with the events
		#self.fig.canvas.mpl_connect('forward_event', self.handle_forward)
		#self.fig.canvas.mpl_connect('backward_event', self.handle_backward)
		self.mpl_toolbar = NavigationToolbar(self.canvas, self.main_frame)
		#self.mpl_toolbar.home_event = self.on_draw
		self.canvas.mpl_connect('home_event', self.on_draw)
		#self.mpl_toolbar.forward = new_forward
		#self.mpl_toolbar.back = new_backward

		# Vertical box layout in the window and setting it up
		vbox = PyQt5.QtWidgets.QVBoxLayout()		
		vbox.addLayout(hlay)
		vbox.addLayout(hlay2)
		vbox.addWidget(self.canvas)  # the matplotlib canvas
		vbox.addWidget(self.mpl_toolbar)
		
		# Adding labels at the bottom
		windowLabel = PyQt5.QtWidgets.QLabel(bottomWindowLabel, self)
		#windowLabel.setOpenExternalLinks(True)
		windowLabel.setAlignment(PyQt5.QtCore.Qt.AlignRight | PyQt5.QtCore.Qt.AlignVCenter)
		vbox.addWidget(windowLabel)

		# We are done...
		self.main_frame.setLayout(vbox)
		self.setCentralWidget(self.main_frame)
		
	"""
	Draws or redraws the plot
	"""
	def on_draw(self, event=None):
		#print ("We are in there")
		# If we want to keep track of the zoom, we save the current zoom
		if (not self.dounzoom):
			left, right = self.axes.get_xlim()
			bottom, top = self.axes.get_ylim()
		# Clear plot
		self.fig.clear()
		self.axes = self.fig.add_subplot(111)
		# Make sure we have data. If not display a message
		if (self.nEta <= 0):
			self.fig.clear()
			self.axes = self.fig.add_subplot(111)
			# Add a label
			self.axes.annotate('Please load data', xy=(.5, 0.5), xycoords='axes fraction', horizontalalignment='center', verticalalignment='center', fontsize=16)
			# Title and labels
			self.axes.set_xlabel("2 theta (degrees)")
			self.axes.set_ylabel("Intensity")
			title = ""
			self.axes.set_title(title, loc='left')
			self.delPointsButton.setDisabled(True)
			# Ready to draw
			self.canvas.draw()
			return
		# Getting plot data
		data = numpy.asarray(self.esgData["data"][self.etaToPlot])
		if (data.size > 0):
			twotheta = data[:,0]
			intensity = data[:,2]
			extralabel = ""
			self.delPointsButton.setDisabled(False)
		else: # empty data for this range
			twotheta = []
			intensity = []
			extralabel = " (no data)"
			self.delPointsButton.setDisabled(True)
			self.axes.annotate('No data for this azimuth', xy=(.5, 0.5), xycoords='axes fraction', horizontalalignment='center', verticalalignment='center', fontsize=12)
		
		# Plot data
		g1 = self.axes.scatter(twotheta, intensity, s=4,  marker='o', facecolors='r', edgecolors='r')
		
		## If enough points, auto background on, and no background yet, we try to generate those points
		if ((len(twotheta)>10*self.nautobg) and self.doautobg and (len(self.xbg) == 0)):
			# Auto-bg points are calculated using the median intensity over a given 2theta range
			# I also-remove the intensity which are too low (below 2/3 of the overall median) to avoid detector dips
			nbg = self.nautobg
			mintt = min(twotheta)
			windowbg = (max(twotheta)-mintt)/nbg
			medI = scipy.median(intensity)
			cutoff = 2.*medI/3.
			#self.xbg.append(twotheta[0])
			#self.ybg.append(intensity[0])
			for i in range(0,nbg):
				x0 = mintt + windowbg*i
				x1 = mintt + windowbg*(i+1)
				ii = []
				test = True
				i = 0
				for tt in twotheta:
					if ((tt>x0) and (tt<x1)):
						if (intensity[i] > cutoff):
							ii.append(intensity[i])
					i += 1
				if (len(ii)>2):
					self.xbg.append((x0+x1)/2.)
					self.ybg.append(scipy.median(ii))
			#self.xbg.append(twotheta[len(twotheta)-1])
			#self.ybg.append(intensity[len(intensity)-1])
			
		
		# If we have background data, add it on
		if (len(self.xbg) > 0):
			g1 = self.axes.plot(self.xbg, self.ybg, color='blue', marker='o', linestyle='solid', linewidth=2, markersize=8)
			#g1 = self.axes.scatter(self.xbg, self.ybg, linestyle='-', s=8,  marker='o', facecolors='b', edgecolors='b')
		
		# Rezoom if need be
		if (not self.dounzoom):
			self.axes.set_xlim(left, right)
			self.axes.set_ylim(bottom, top)
			self.dounzoom = True
		
		# Title and labels
		self.axes.set_xlabel("2theta")
		self.axes.set_ylabel("intensity")
		title = "Id %d, eta %s %s" % (self.etaToPlot, self.esgData["etas"][self.etaToPlot], extralabel)
		self.axes.set_title(title, loc='left')
		
		# Ready to draw
		self.canvas.draw()

	"""
	Event processing: we need to change dataset based on text input
	"""
	def new_eta(self):
		if (self.nEta>0):
			try:
				i = int(self.etaNBox.text())
				self.etaToPlot = (i) % self.nEta
				self.etaNBox.setText("%d" % (self.etaToPlot))
				self.xbg = []				# Used for creating background
				self.ybg = []				# Used for creating background
				self.subtractBgButton.setDisabled(True)
				self.on_draw()
			except ValueError:
				PyQt5.QtWidgets.QMessageBox.critical(self, "Error", "Not an integer")


	"""
	Event processing when left arrow is click (move to previous dataset)
	"""
	def handle_backward(self,evt=None):
		if (self.nEta>0):
			self.etaToPlot = (self.etaToPlot-1) % self.nEta
			self.etaNBox.setText("%d" % (self.etaToPlot))
			self.xbg = []				# Used for creating background
			self.ybg = []				# Used for creating background
			self.subtractBgButton.setDisabled(True)
			self.on_draw()

	"""
	Event processing when right arrow is click (move to next dataset)
	"""
	def handle_forward(self,evt=None):
		if (self.nEta>0):
			self.etaToPlot = (self.etaToPlot+1) % self.nEta
			self.etaNBox.setText("%d" % (self.etaToPlot))
			self.xbg = []				# Used for creating background
			self.ybg = []				# Used for creating background
			self.subtractBgButton.setDisabled(True)
			self.on_draw()
		
	"""
	Event processing to remove data points from a dataset
	"""
	def remove_points(self,evt=None):
		if (self.nEta > 0):
			# Getting the X and Y ranges to be remove
			left, right = self.axes.get_xlim()
			bottom, top = self.axes.get_ylim()
			data = numpy.asarray(self.esgData["data"][self.etaToPlot])
			# Saving old data for undos
			self.olddata.append(data)
			self.oldetaToPlot.append(self.etaToPlot)
			# Saving range to remove to mask
			self.mask.append({"set":True, "eta": self.etaToPlot, "clear2thetamin": left, "clear2thetamax": right})
			# Remove points within our range
			toremove = []
			for i in range(data.shape[0]):
				twotheta = data[i,0]
				intensity = data[i,2]
				if ((twotheta>left) and (twotheta<right) and (intensity>bottom) and (intensity<top)):
					toremove.append(i)
			data = numpy.delete(data,toremove,0)
			self.esgData["data"][self.etaToPlot] = data.tolist()
			# We change something. We should ask confirmation for saving before closing the app, we can also undo stuff from now on
			self.needToSave = True
			self.cancelButton.setDisabled(False)
			# Delete the background to avoid strange effects
			self.xbg = []				# Used for creating background
			self.ybg = []				# Used for creating background
			# Redraw everything
			self.on_draw()
			
	"""
	Restrict 2 theta range based on user input
	"""
	def edit_twothetarange(self,evt=None):
		if (self.nEta > 0):
			test = TThethaRangeDialog(self)
			result = test.exec_()
			if (test.isOk()):
				min2theta,max2theta = test.getInputs()
				# Saving old data for undos
				data = copy.deepcopy(self.esgData["data"])
				self.olddata.append(data)
				self.oldetaToPlot.append("all")
				# Remove points outside our 2theta range
				for i in range(0,self.nEta):
					thisetadata = numpy.asarray(self.esgData["data"][i])
					toremove = []
					for j in range(thisetadata.shape[0]):
						twotheta = thisetadata[j,0]
						intensity = thisetadata[j,2]
						if ((twotheta>max2theta) or (twotheta<min2theta)):
							toremove.append(j)
					thisetadata = numpy.delete(thisetadata,toremove,0)
					self.esgData["data"][i] = thisetadata.tolist()
				self.needToSave = True
				self.cancelButton.setDisabled(False)
				self.on_draw()
		return
	
	"""
	Event processing when we want to cancel. Go back to the last version.
	"""
	def cancel_last(self,evt=None):
		if ((len(self.olddata)>0) and (self.nEta > 0)):
			# We pop the elements of the self.oldetaToPlot and self.olddata lists and set them as the new data
			test = self.oldetaToPlot.pop()
			if (test == "all"): # We changed the 2 theta range or something that affects all azimuths
				self.esgData["data"] = self.olddata.pop()
			else: # We changed something that affects a single azimuth
				self.etaToPlot = test
				self.etaNBox.setText("%d" % (self.etaToPlot))
				self.esgData["data"][self.etaToPlot] = self.olddata.pop()
				# Replot
			if (len(self.olddata) == 0):
				self.cancelButton.setDisabled(True)
			self.mask.pop() # Remove last addition to mask
			self.on_draw()
	
	"""
	Event to quit the app
	"""
	def closeEvent(self,evt=None):
		if (self.needToSave):
			buttonReply = PyQt5.QtWidgets.QMessageBox.question(self, 'Data not saved', "Data not saved. Quit anyway?", PyQt5.QtWidgets.QMessageBox.Yes | PyQt5.QtWidgets.QMessageBox.No, PyQt5.QtWidgets.QMessageBox.No)
			if (buttonReply == PyQt5.QtWidgets.QMessageBox.No):
				if (isinstance(evt,PyQt5.QtGui.QCloseEvent)):
					evt.ignore()
				return
		if (isinstance(evt,PyQt5.QtGui.QCloseEvent)):
			evt.accept()
		sys.exit(2)

	"""
	Save current dataset to ESG format
	"""
	def save_esg(self,evt=None):
		options = PyQt5.QtWidgets.QFileDialog.Options()
		fileName, _ = PyQt5.QtWidgets.QFileDialog.getSaveFileName(self,"Save new data as...", "","Esg Files (*.esg);;All Files (*)", options=options)
		if fileName:
			saveEsgToFile(self.esgData,fileName)
			self.needToSave = False
			path, name = os.path.split(fileName)
			self.title = "MAUD ESG edit: " + name
			self.setWindowTitle(self.title)

	"""
	Save current mask
	"""
	def save_mask(self,evt=None):
		options = PyQt5.QtWidgets.QFileDialog.Options()
		fileName, _ = PyQt5.QtWidgets.QFileDialog.getSaveFileName(self,"Save mask as...", "","maudESGEdit Mask Files (*.msk);;All Files (*)", options=options)
		if fileName:
			saveMaskToFile(self.mask,fileName)
			
	"""
	Load and apply a mask
	"""
	def load_and_apply_mask(self,evt=None):
		if (self.nEta<1):
			PyQt5.QtWidgets.QMessageBox.critical(self, "Error", "Please load some data first")
			return
		if (self.needToSave):
			buttonReply = PyQt5.QtWidgets.QMessageBox.question(self, 'Data not saved', "Data not saved. Load and apply a mask anyway? Previous mask will be lost. Things may crash if you try to undo.", PyQt5.QtWidgets.QMessageBox.Yes | PyQt5.QtWidgets.QMessageBox.No, PyQt5.QtWidgets.QMessageBox.No)
			if (buttonReply == PyQt5.QtWidgets.QMessageBox.No):
				return
		options = PyQt5.QtWidgets.QFileDialog.Options()
		filename, _ = PyQt5.QtWidgets.QFileDialog.getOpenFileName(self,"Select an mask file...", self.pathtomask,"maudESGEdit Mask Files (*.msk);;All Files (*)", options=options)
		if filename:
			self.pathtomask = os.path.dirname(filename)
			self.mask = loadMaskFromFile(filename)
			# Saving old data for undos
			data = copy.deepcopy(self.esgData["data"])
			self.olddata.append(data)
			self.oldetaToPlot.append("all")
			# Remove points in the mask
			for item in self.mask:
				if (item["set"]):
					#self.mask.append({"set":True, "eta": self.etaToPlot, "clear2thetamin": left, "clear2thetamax": right})
					eta = item["eta"]
					min2theta = item["clear2thetamin"]
					max2theta = item["clear2thetamax"]
					thisetadata = numpy.asarray(self.esgData["data"][eta])
					toremove = []
					for j in range(thisetadata.shape[0]):
						twotheta = thisetadata[j,0]
						intensity = thisetadata[j,2]
						if ((twotheta<max2theta) and (twotheta>min2theta)):
							toremove.append(j)
					thisetadata = numpy.delete(thisetadata,toremove,0)
					self.esgData["data"][eta] = thisetadata.tolist()
			self.needToSave = True
			self.cancelButton.setDisabled(False)
			self.on_draw()
			
			
	"""
	Open a different esg
	"""
	def open_esg(self,evt=None):
		if (self.needToSave):
			buttonReply = PyQt5.QtWidgets.QMessageBox.question(self, 'Data not saved', "Data not saved. Load a new dataset anyway?", PyQt5.QtWidgets.QMessageBox.Yes | PyQt5.QtWidgets.QMessageBox.No, PyQt5.QtWidgets.QMessageBox.No)
			if (buttonReply == PyQt5.QtWidgets.QMessageBox.No):
				return
		options = PyQt5.QtWidgets.QFileDialog.Options()
		filename, _ = PyQt5.QtWidgets.QFileDialog.getOpenFileName(self,"Select an ESG file...", "","Esg Files (*.esg);;All Files (*)", options=options)
		if filename:
			self.esgData = parseESG(filename)
			path, name = os.path.split(filename)
			self.title = "MAUD ESG edit: " + name
			self.filename = filename
			self.fileSaveHint = filename
			self.nEta = len(self.esgData["etas"])
			self.etaLabel.setText("Spectrum (0-%d) : " % (self.nEta-1))
			self.etaToPlot = 0
			self.etaNBox.setText("%d" % (self.etaToPlot))
			self.olddata = [] # Deleting cached old data to avoid confusion
			self.mask = []	# clear mask
			self.setWindowTitle(self.title)
			self.needToSave = False
			self.tthetaButton.setDisabled(False)
			self.subtractBgButton.setDisabled(True)
			self.xbg = []				# Used for creating background
			self.ybg = []				# Used for creating background
			self.on_draw()
		#else:
		#	PyQt5.QtWidgets.QMessageBox.critical(self, "Error", "File opening failed")
	
	"""
	Opens the about window
	"""
	def about(self,evt=None):
		dialog = textWindow("About MAUD ESG editor...", aboutWindowText, self)
		dialog.exec_()
		return
	
	"""
	Opens the help window
	"""
	def helpwindow(self,evt=None):
		dialog = textWindow("Help for MAUD ESG editor...", helpWindowText, self)
		dialog.exec_()
		return
	
	"""
	The user clicked somewhere, if right button, we generate background points
	"""
	def on_press(self, event):
		if (int(event.button) == 3):
			self.xbg.append(event.xdata)
			self.ybg.append(event.ydata)
			# Sort by increasing x-values
			#yx = zip(self.xbg, self.ybg)
			#yx.sort()
			#self.xbg = [x for x, y in yx]
			#self.ybg = [y for x, y in yx]
			if (len(self.xbg)>1):
				# Sorting background points by increasing 2theta
				xcopy = self.xbg[:]
				ycopy = self.ybg[:]
				indices = sorted(range(len(xcopy)), key=lambda k: xcopy[k])
				self.xbg = []
				self.ybg = []
				for i in indices:
					self.xbg.append(xcopy[i])
					self.ybg.append(ycopy[i])
				# Allowing for background subtraction
				self.subtractBgButton.setDisabled(False)
			# Ready to plot
			self.dounzoom = False
			self.on_draw()
	
	"""
	Delete all background points
	"""
	def remove_bg_points(self, event=None):
		self.xbg = []
		self.ybg = []
		self.subtractBgButton.setDisabled(True)
		self.on_draw()
	
	"""
	We need to subtract the background from the data
	"""
	def subtract_background(self, event):
		if (len(self.xbg)>0):
			bg = interpolate.interp1d(self.xbg, self.ybg)
			data = numpy.asarray(self.esgData["data"][self.etaToPlot])
			# Saving old data for undos
			self.olddata.append(data.copy())
			self.oldetaToPlot.append(self.etaToPlot)
			self.mask.append({"set":False}); # Adding an empty value in saved mask. Necessary for proper handle of undos
			# Remove background within our range
			xmin = min(self.xbg)
			xmax = max(self.xbg)
			for i in range(data.shape[0]):
				twotheta = data[i,0]
				intensity = data[i,2]
				if ((twotheta>xmin) and (twotheta<xmax)):
					try:
						data[i,2] = data[i,2] - bg(twotheta)
					except TypeError:
						print ("Small error in background subtraction, we keep going")
			self.esgData["data"][self.etaToPlot] = data.tolist()
			# We change something. We should ask confirmation for saving before closing the app, we can also undo stuff from now on
			self.needToSave = True
			self.cancelButton.setDisabled(False)
			self.xbg=[]
			self.ybg=[]
			self.subtractBgButton.setDisabled(True)
			# Redraw everything
			self.on_draw()
			
	"""
	Add a constant value for all intensities at a given azimuth
	"""
	def shift_data(self, event):
		if (self.nEta <= 0):
			return
		shift,ok = PyQt5.QtWidgets.QInputDialog.getDouble(self,"Shift data by","How much shall we add to intensities (current azimuth only)")
		if ok:
			data = numpy.asarray(self.esgData["data"][self.etaToPlot])
			if (data.size > 0):
				# Saving old data for undos
				self.olddata.append(data.copy())
				self.oldetaToPlot.append(self.etaToPlot)
				self.mask.append({"set":False}); # Adding an empty value in saved mask. Necessary for proper handle of undos
				# Adding to the intensity
				data[:,2] += shift
				self.esgData["data"][self.etaToPlot] = data.tolist()
				# We change something. We should ask confirmation for saving before closing the app, we can also undo stuff from now on
				self.needToSave = True
				self.cancelButton.setDisabled(False)
				self.xbg=[]
				self.ybg=[]
				# Redraw everything
				self.on_draw()

	"""
	Add a constant value for all intensities at all azimuth
	"""
	def shift_data_all(self, event):
		if (self.nEta <= 0):
			return
		shift,ok = PyQt5.QtWidgets.QInputDialog.getDouble(self,"Shift data by","How much shall we add to intensities (all azimuthal angles)")
		if ok:
			# Saving old data for undos
			data = copy.deepcopy(self.esgData["data"])
			self.olddata.append(data)
			self.oldetaToPlot.append("all")
			self.mask.append({"set":False}); # Adding an empty value in saved mask. Necessary for proper handle of undos
			# Shifting intensities
			for i in range(0,self.nEta):
				thisetadata = numpy.asarray(self.esgData["data"][i])
				if (thisetadata.size > 0):
					# Adding to the intensity
					thisetadata[:,2] += shift
					self.esgData["data"][i] = thisetadata.tolist()
			# We change something. We should ask confirmation for saving before closing the app, we can also undo stuff from now on
			self.needToSave = True
			self.cancelButton.setDisabled(False)
			self.xbg=[]
			self.ybg=[]
			# Redraw everything
			self.on_draw()

	"""
	Set a fixed minimum intensity for all azimuth
	"""
	def setmin_data_all(self, event):
		if (self.nEta <= 0):
			return
		minval,ok = PyQt5.QtWidgets.QInputDialog.getDouble(self,"Minimum","Minimum intensity to set at all azimuthal angles")
		if ok:
			# Saving old data for undos
			data = copy.deepcopy(self.esgData["data"])
			self.olddata.append(data)
			self.oldetaToPlot.append("all")
			self.mask.append({"set":False}); # Adding an empty value in saved mask. Necessary for proper handle of undos
			# Shifting intensities
			for i in range(0,self.nEta):
				thisetadata = numpy.asarray(self.esgData["data"][i])
				if (thisetadata.size > 0):
					# Adding to the intensity
					shift = minval-min(thisetadata[:,2])
					thisetadata[:,2] += shift
					self.esgData["data"][i] = thisetadata.tolist()
			# We change something. We should ask confirmation for saving before closing the app, we can also undo stuff from now on
			self.needToSave = True
			self.cancelButton.setDisabled(False)
			self.xbg=[]
			self.ybg=[]
			# Redraw everything
			self.on_draw()
	
	"""
	Turn on or off the auto-background feature
	"""
	def changeautobg(self, state):
		if (state == PyQt5.QtCore.Qt.Checked):
			self.doautobg = True
			self.autobgNBox.setDisabled(False)
		else:
			self.doautobg = False
			self.autobgNBox.setDisabled(True)
		self.xbg=[]
		self.ybg=[]
		# Redraw everything
		self.on_draw()
	
	"""
	Change the number of autobackground points
	"""
	def changenautobg(self):
		try:
			self.nautobg = int(self.autobgNBox.text())
			self.xbg=[]
			self.ybg=[]
			# Redraw everything
			self.on_draw()
		except ValueError:
			PyQt5.QtWidgets.QMessageBox.critical(self, "Error", "Not an integer")

#################################################################
#
# Application icon in XPM format, encoded in base64
#
#################################################################

iconXPMbase64 = "LyogWFBNICovCnN0YXRpYyBjaGFyICogbWF1ZEVTR0VkaXRfeHBtW10gPSB7CiI5NSA4NiAxOTM5IDIiLAoiICAJYyBOb25lIiwKIi4gCWMgIzAwMDBDRiIsCiIrIAljICMwMDAwQ0UiLAoiQCAJYyAjMUExQUQzIiwKIiMgCWMgIzVDNUNERiIsCiIkIAljICNBM0EzRUUiLAoiJSAJYyAjQzZDNkY0IiwKIiYgCWMgI0U0RTRGOSIsCiIqIAljICNFN0U3RkIiLAoiPSAJYyAjQ0VDRUY1IiwKIi0gCWMgI0FFQUVFRiIsCiI7IAljICM2RTZFRTMiLAoiPiAJYyAjMjkyOUQ2IiwKIiwgCWMgIzAxMDFDRSIsCiInIAljICMwMDAwQ0QiLAoiKSAJYyAjMUQxREQ0IiwKIiEgCWMgIzk5OTlFQyIsCiJ+IAljICNGRkZGRkYiLAoieyAJYyAjQkNCQ0YzIiwKIl0gCWMgIzMyMzJEOCIsCiJeIAljICM1NDU0REUiLAoiLyAJYyAjRjlGOUZFIiwKIiggCWMgIzgyODJFNyIsCiJfIAljICMwNDA0Q0YiLAoiOiAJYyAjNzc3N0U1IiwKIjwgCWMgI0FDQUNFRiIsCiJbIAljICMwNzA3Q0YiLAoifSAJYyAjNUE1QURGIiwKInwgCWMgIzk1OTVFQSIsCiIxIAljICMyNzI3RDUiLAoiMiAJYyAjRkJGQkZFIiwKIjMgCWMgI0FBQUFFRiIsCiI0IAljICNGRkZFRkUiLAoiNSAJYyAjRkVGQkZCIiwKIjYgCWMgI0Y4RjRGNCIsCiI3IAljICNGMUVGRUQiLAoiOCAJYyAjRTlFOUU2IiwKIjkgCWMgI0U0RTZFMyIsCiIwIAljICNFMUU1RTAiLAoiYSAJYyAjREVFMURCIiwKImIgCWMgI0Q5RERENyIsCiJjIAljICNEN0Q5RDUiLAoiZCAJYyAjRDZEOEQ0IiwKImUgCWMgI0Q3RDlENiIsCiJmIAljICNEQURCRDgiLAoiZyAJYyAjREVERkRFIiwKImggCWMgI0UyRTJFMSIsCiJpIAljICNFNkU0RTMiLAoiaiAJYyAjRUZFQkVBIiwKImsgCWMgI0Y3RjBGMCIsCiJsIAljICNGQ0Y1RjciLAoibSAJYyAjRkNGQUZBIiwKIm4gCWMgI0ZFRkVGRSIsCiJvIAljICNFMkUyRjkiLAoicCAJYyAjMEEwQUQwIiwKInEgCWMgIzJCMkJENiIsCiJyIAljICNGRUZDRkMiLAoicyAJYyAjRjRFRUVFIiwKInQgCWMgI0VBRENEQyIsCiJ1IAljICNFMUNGRDAiLAoidiAJYyAjRDJDMUJGIiwKIncgCWMgI0JFQjZCMyIsCiJ4IAljICNCQkJEQjMiLAoieSAJYyAjQkNDNUI3IiwKInogCWMgI0JBQzhCNiIsCiJBIAljICNCRENDQjkiLAoiQiAJYyAjQkFDN0I2IiwKIkMgCWMgI0I2QzNCMiIsCiJEIAljICNCNEJGQjEiLAoiRSAJYyAjQjJCQ0IwIiwKIkYgCWMgI0IwQjdBRSIsCiJHIAljICNBRkI0QUUiLAoiSCAJYyAjQjFCM0FGIiwKIkkgCWMgI0I0QjBBRCIsCiJKIAljICNCRkIxQjAiLAoiSyAJYyAjRDVCOEI5IiwKIkwgCWMgI0UzQzRDNiIsCiJNIAljICNFNkQwRDAiLAoiTiAJYyAjRUZFMUUxIiwKIk8gCWMgI0Y2RURFRCIsCiJQIAljICNGOUYyRjMiLAoiUSAJYyAjNUQ1REUwIiwKIlIgCWMgIzcxNzFFNCIsCiJTIAljICNGNkY1RjUiLAoiVCAJYyAjRTVFMEUwIiwKIlUgCWMgI0REQ0ZEMCIsCiJWIAljICNEOUMwQzMiLAoiVyAJYyAjRDRCNUI4IiwKIlggCWMgI0Q0QjFCMiIsCiJZIAljICNEMkFFQUUiLAoiWiAJYyAjQzhBREFDIiwKImAgCWMgI0I2QUJBQiIsCiIgLgljICNCNUI1QjAiLAoiLi4JYyAjQkRDNkI5IiwKIisuCWMgI0JGQ0VCQiIsCiJALgljICNDM0QyQkMiLAoiIy4JYyAjQzFDRUI5IiwKIiQuCWMgI0JEQ0FCNyIsCiIlLgljICNCQUM3QjciLAoiJi4JYyAjQjhDMkI1IiwKIiouCWMgI0I0QkFCMiIsCiI9LgljICNBRkFGQUQiLAoiLS4JYyAjQjBBQ0FCIiwKIjsuCWMgI0JDQUJBQiIsCiI+LgljICNEM0FEQUUiLAoiLC4JYyAjRTBBREFFIiwKIicuCWMgI0UxQjBBRiIsCiIpLgljICNFMUI0QjMiLAoiIS4JYyAjRTJCQkJBIiwKIn4uCWMgI0VFRDFEMCIsCiJ7LgljICNGMUUwREUiLAoiXS4JYyAjRkFGM0YyIiwKIl4uCWMgI0IzQjNGMCIsCiIvLgljICNDMkMyRjQiLAoiKC4JYyAjRURFN0U3IiwKIl8uCWMgI0U0RDZENiIsCiI6LgljICNEM0I4QjkiLAoiPC4JYyAjRDJCMEIxIiwKIlsuCWMgI0Q1QUNBRSIsCiJ9LgljICNEOEFCQUQiLAoifC4JYyAjRDlBQkFDIiwKIjEuCWMgI0RBQUJBQiIsCiIyLgljICNENkFEQUQiLAoiMy4JYyAjQ0JBREFEIiwKIjQuCWMgI0JGQUNBQiIsCiI1LgljICNCOEFDQUIiLAoiNi4JYyAjQjZCMUFCIiwKIjcuCWMgI0IxQUVBQyIsCiI4LgljICM3MjcwQTgiLAoiOS4JYyAjM0QzQ0E1IiwKIjAuCWMgIzE4MTdBMyIsCiJhLgljICMwODA4QTEiLAoiYi4JYyAjMEIwQUEyIiwKImMuCWMgIzI0MjFBMyIsCiJkLgljICM0MzNDQTUiLAoiZS4JYyAjODQ3M0E5IiwKImYuCWMgI0NDQTZBRCIsCiJnLgljICM2OTUwQTciLAoiaC4JYyAjMzAyNEE0IiwKImkuCWMgIzhGNjlBOCIsCiJqLgljICNFOUFDQUIiLAoiay4JYyAjRTdBREFEIiwKImwuCWMgI0UzQUZBRSIsCiJtLgljICNFMkI4QjYiLAoibi4JYyAjRUFDOUM3IiwKIm8uCWMgI0Y2RTJFMCIsCiJwLgljICNGQkY0RjQiLAoicS4JYyAjQkJCQkU2IiwKInIuCWMgIzYxNjFDNSIsCiJzLgljICMyOTI5QjAiLAoidC4JYyAjMTQxNEE4IiwKInUuCWMgIzA1MDVBMyIsCiJ2LgljICM0RTRFQkUiLAoidy4JYyAjQTJBMkREIiwKInguCWMgI0Y4RjhGQyIsCiJ5LgljICM4Nzg3RDMiLAoiei4JYyAjM0EzQUI3IiwKIkEuCWMgIzlDOUNEQiIsCiJCLgljICNGQ0ZDRkYiLAoiQy4JYyAjMEIwQkQwIiwKIkQuCWMgI0VFRUVGQyIsCiJFLgljICM4RjhGRDYiLAoiRi4JYyAjMTkxOUFBIiwKIkcuCWMgIzAwMDBBMSIsCiJILgljICMwMjAxQTEiLAoiSS4JYyAjNDkzQUE1IiwKIkouCWMgI0RDQUJBQiIsCiJLLgljICNEREFCQUMiLAoiTC4JYyAjREZBQkFDIiwKIk0uCWMgI0Q4QURBRCIsCiJOLgljICNENEFFQUUiLAoiTy4JYyAjODk3MkFBIiwKIlAuCWMgIzEyMEZBMiIsCiJRLgljICMwOTA3QTIiLAoiUi4JYyAjMjkxRUEzIiwKIlMuCWMgI0VEQUJBQyIsCiJULgljICNFRUFCQUMiLAoiVS4JYyAjRURBQ0FDIiwKIlYuCWMgI0VCQURBRSIsCiJXLgljICNFNkI1QjYiLAoiWC4JYyAjRTlDQ0NCIiwKIlkuCWMgI0Y5RjBFRiIsCiJaLgljICNDNkM2RUEiLAoiYC4JYyAjMzMzM0IzIiwKIiArCWMgIzBCMEJBNSIsCiIuKwljICMwMjAyQTIiLAoiKysJYyAjMkUyRUIyIiwKIkArCWMgIzIyMjJENCIsCiIjKwljICMwOTA5Q0YiLAoiJCsJYyAjNDU0NUJCIiwKIiUrCWMgI0MzOUFBQSIsCiImKwljICNEQUFDQUIiLAoiKisJYyAjRENBQ0FCIiwKIj0rCWMgI0UxQUJBQyIsCiItKwljICNFMEFDQUMiLAoiOysJYyAjNjk1MUE2IiwKIj4rCWMgIzA0MDNBMSIsCiIsKwljICMyMjFBQTIiLAoiJysJYyAjNDQzM0E0IiwKIikrCWMgIzQ1MzNBNCIsCiIhKwljICMxMDBCQTEiLAoifisJYyAjRjBBQkFDIiwKInsrCWMgI0YxQUJBQyIsCiJdKwljICNGMkFCQUMiLAoiXisJYyAjRjFBQ0FDIiwKIi8rCWMgI0VFQUNBQyIsCiIoKwljICNFN0IwQjEiLAoiXysJYyAjRTNDN0M3IiwKIjorCWMgI0YwRThFOCIsCiI8KwljICM4RDhERDUiLAoiWysJYyAjMDMwM0EyIiwKIn0rCWMgIzBEMERBNiIsCiJ8KwljICMyRjJGQjMiLAoiMSsJYyAjM0YzRkI4IiwKIjIrCWMgIzA5MDlBNSIsCiIzKwljICMxODE4QUEiLAoiNCsJYyAjMzczN0Q4IiwKIjUrCWMgIzBGMEZEMSIsCiI2KwljICNBRkFGRTIiLAoiNysJYyAjNDE0MUI5IiwKIjgrCWMgIzI4MjhBRiIsCiI5KwljICMxRjFGQUMiLAoiMCsJYyAjMjIyMkFFIiwKImErCWMgIzI0MjRBRCIsCiJiKwljICMyOTI2QTkiLAoiYysJYyAjMjkyM0E2IiwKImQrCWMgIzI5MjFBNCIsCiJlKwljICMyOTIwQTMiLAoiZisJYyAjMTAwQ0EyIiwKImcrCWMgI0JCOTFBQSIsCiJoKwljICNEQkFDQUIiLAoiaSsJYyAjRENBQ0FDIiwKImorCWMgI0RFQUJBQiIsCiJrKwljICNFMEFCQUIiLAoibCsJYyAjRTJBQkFCIiwKIm0rCWMgIzhDNkFBNyIsCiJuKwljICM1NzQwQTQiLAoibysJYyAjQ0Q5OEFBIiwKInArCWMgI0U4QUJBQyIsCiJxKwljICNFQUFCQUMiLAoicisJYyAjRUJBQkFDIiwKInMrCWMgI0MyOEVBQSIsCiJ0KwljICMzRTJEQTQiLAoidSsJYyAjRUZBQkFEIiwKInYrCWMgI0YyQUNBRCIsCiJ3KwljICNGM0FDQUMiLAoieCsJYyAjRjRBQ0FDIiwKInkrCWMgI0VGQURBRSIsCiJ6KwljICNEOUIwQjAiLAoiQSsJYyAjNzQ2REFFIiwKIkIrCWMgIzAxMDFBMSIsCiJDKwljICMxRDFEQUIiLAoiRCsJYyAjOUQ5RERCIiwKIkUrCWMgI0VDRUNGOCIsCiJGKwljICM3MTcxQ0IiLAoiRysJYyAjM0YzRkRBIiwKIkgrCWMgI0RFREVGNCIsCiJJKwljICNGM0VERUQiLAoiSisJYyAjREVDOUM5IiwKIksrCWMgI0Q3QjJCMyIsCiJMKwljICNEQ0FEQUUiLAoiTSsJYyAjREVBQ0FCIiwKIk4rCWMgI0RGQUNBQiIsCiJPKwljICM1NzQzQTUiLAoiUCsJYyAjQkM5MUFCIiwKIlErCWMgI0RFQUJBQyIsCiJSKwljICNFMUFCQUIiLAoiUysJYyAjRTNBQkFCIiwKIlQrCWMgI0REQTZBQiIsCiJVKwljICMxMTBEQTIiLAoiVisJYyAjNjg0REE1IiwKIlcrCWMgI0U3QUJBQiIsCiJYKwljICNFOEFCQUIiLAoiWSsJYyAjRTlBQkFDIiwKIlorCWMgI0VDQUJBQyIsCiJgKwljICNFOUFBQUIiLAoiIEAJYyAjNDczM0E0IiwKIi5ACWMgI0UzQTJBQyIsCiIrQAljICNGM0FFQUQiLAoiQEAJYyAjRjVBRUFFIiwKIiNACWMgI0Y2QUVBRiIsCiIkQAljICNGNkFEQUYiLAoiJUAJYyAjRjdBREFEIiwKIiZACWMgI0Y0QUVBQyIsCiIqQAljICNCODhDQUIiLAoiPUAJYyAjMDQwNEExIiwKIi1ACWMgIzNDM0NCNyIsCiI7QAljICNGNkY2RkMiLAoiPkAJYyAjQTZBNkRFIiwKIixACWMgIzA4MDhBNCIsCiInQAljICNGQ0ZDRkMiLAoiKUAJYyAjRjBFQUVBIiwKIiFACWMgI0RFQzJDNCIsCiJ+QAljICNEOUFGQjIiLAoie0AJYyAjREVBQ0FFIiwKIl1ACWMgI0UwQUNBQiIsCiJeQAljICM1ODQzQTUiLAoiL0AJYyAjQkU5MUFBIiwKIihACWMgI0UwQUJBQyIsCiJfQAljICNFNEFCQUIiLAoiOkAJYyAjRTVBQkFCIiwKIjxACWMgIzlDNzVBOCIsCiJbQAljICMxNDBGQTIiLAoifUAJYyAjRTNBOEFCIiwKInxACWMgI0U4QUNBQiIsCiIxQAljICNFOUFCQUIiLAoiMkAJYyAjRUFBQkFCIiwKIjNACWMgI0VCQUJBQiIsCiI0QAljICNFREFCQUIiLAoiNUAJYyAjRUZBQ0FDIiwKIjZACWMgI0NBOTFBQSIsCiI3QAljICNEMzk4QUIiLAoiOEAJYyAjRjZBRkFFIiwKIjlACWMgI0Y4QUZBRiIsCiIwQAljICNGOUFGQjAiLAoiYUAJYyAjRkFCMEIyIiwKImJACWMgI0ZBQjBCMCIsCiJjQAljICMyNTFCQTMiLAoiZEAJYyAjMjMyMkE0IiwKImVACWMgI0MyQzJDNiIsCiJmQAljICNGMEYwRjAiLAoiZ0AJYyAjNUQ1REMzIiwKImhACWMgI0VCRTJFMiIsCiJpQAljICNEQ0JGQkYiLAoiakAJYyAjRERBREFGIiwKImtACWMgI0UxQUNBRSIsCiJsQAljICNFMkFDQUMiLAoibUAJYyAjRTFBQ0FCIiwKIm5ACWMgIzc1NTdBNyIsCiJvQAljICM0NjM0QTQiLAoicEAJYyAjRUJBQkFEIiwKInFACWMgI0VDQUJBQiIsCiJyQAljICNFRkFCQUMiLAoic0AJYyAjRUZBQ0FEIiwKInRACWMgI0YwQUNBRSIsCiJ1QAljICM0OTM0QTUiLAoidkAJYyAjMEQwQUEyIiwKIndACWMgI0VFQUJCMCIsCiJ4QAljICNGQUIxQjEiLAoieUAJYyAjRkJCMkIyIiwKInpACWMgI0ZEQjNCNCIsCiJBQAljICNGREI0QjQiLAoiQkAJYyAjQTY3NUFFIiwKIkNACWMgIzhCOEFBQiIsCiJEQAljICNBQkFDQUQiLAoiRUAJYyAjQkZCRUJFIiwKIkZACWMgI0YxRUZFRiIsCiJHQAljICNFOUU5RjYiLAoiSEAJYyAjMUMxQ0FDIiwKIklACWMgI0U2RENEQyIsCiJKQAljICNEREMxQzAiLAoiS0AJYyAjREVBRkFFIiwKIkxACWMgI0UzQUNBQyIsCiJNQAljICNDNzk4QUIiLAoiTkAJYyAjOUI3NkE4IiwKIk9ACWMgI0RGQThBQiIsCiJQQAljICM2MjRBQTUiLAoiUUAJYyAjQzg5QUFBIiwKIlJACWMgI0UxQUNBQyIsCiJTQAljICNFM0FCQUMiLAoiVEAJYyAjRTRBQkFDIiwKIlVACWMgIzZCNTBBNiIsCiJWQAljICM0ODM1QTQiLAoiV0AJYyAjRUNBQ0FCIiwKIlhACWMgI0VFQUNBQiIsCiJZQAljICNGMEFEQUQiLAoiWkAJYyAjRjFBREFFIiwKImBACWMgI0UxQTFBRSIsCiIgIwljICM5MjY4QTkiLAoiLiMJYyAjQ0E5MUFFIiwKIisjCWMgI0Y5QjNCMiIsCiJAIwljICNGQ0I0QjQiLAoiIyMJYyAjRkRCNkI1IiwKIiQjCWMgI0ZGQjdCNyIsCiIlIwljICNGRkI4QjgiLAoiJiMJYyAjRkZCOUI4IiwKIiojCWMgIzM2MjdBNiIsCiI9IwljICMyRDJCQTQiLAoiLSMJYyAjQUZBRUFGIiwKIjsjCWMgI0FDQUNBRCIsCiI+IwljICNBREFDQUMiLAoiLCMJYyAjQzJDMkMxIiwKIicjCWMgI0VFRUVFRSIsCiIpIwljICM4QThBRDQiLAoiISMJYyAjMzAzMEIyIiwKIn4jCWMgIzg0ODREMiIsCiJ7IwljICNFQUUyRTMiLAoiXSMJYyAjRERCRUJGIiwKIl4jCWMgI0UwQUVBRCIsCiIvIwljICNFNUFDQUIiLAoiKCMJYyAjRTZBQ0FDIiwKIl8jCWMgIzNBMkNBNCIsCiI6IwljICM4QzY5QTkiLAoiPCMJYyAjMUMxNUEyIiwKIlsjCWMgI0UyQUJBQyIsCiJ9IwljICNFNEFCQUQiLAoifCMJYyAjRTVBQkFEIiwKIjEjCWMgI0U0QUNBQyIsCiIyIwljICM4MzYxQTYiLAoiMyMJYyAjRDQ5QkFCIiwKIjQjCWMgI0VEQUNBQiIsCiI1IwljICNGMkFFQUYiLAoiNiMJYyAjRjVBRkIwIiwKIjcjCWMgI0Y3QjFCMSIsCiI4IwljICNGQ0I1QjQiLAoiOSMJYyAjRkRCN0I3IiwKIjAjCWMgI0ZFQjlCOSIsCiJhIwljICNGRkJCQkIiLAoiYiMJYyAjRkZCQ0JDIiwKImMjCWMgI0VCQUZCQyIsCiJkIwljICM3MDZCQUQiLAoiZSMJYyAjQjZCNEI0IiwKImYjCWMgI0IxQUZCMCIsCiJnIwljICNBRUFEQUQiLAoiaCMJYyAjQUJBQ0FCIiwKImkjCWMgI0MwQzFDMSIsCiJqIwljICNGQUZBRkEiLAoiayMJYyAjREVDOUNBIiwKImwjCWMgI0RBQjJCMiIsCiJtIwljICNFMEFEQUQiLAoibiMJYyAjRTJBREFDIiwKIm8jCWMgI0U1QUNBQyIsCiJwIwljICMxNTBGQTEiLAoicSMJYyAjNkQ1MUE3IiwKInIjCWMgI0NDOUFBQSIsCiJzIwljICNEOUE0QUEiLAoidCMJYyAjRTJBQ0FCIiwKInUjCWMgI0U2QUJBQyIsCiJ2IwljICNFNEFDQUQiLAoidyMJYyAjQkY5MEFCIiwKIngjCWMgIzJGMjJBMyIsCiJ5IwljICNCRjhCQTkiLAoieiMJYyAjRUZBREFDIiwKIkEjCWMgI0YzQUVBRiIsCiJCIwljICNGNUIwQjEiLAoiQyMJYyAjRjhCMkIyIiwKIkQjCWMgI0ZCQjRCNCIsCiJFIwljICNGRkJFQkUiLAoiRiMJYyAjRkZDMUMxIiwKIkcjCWMgI0FBODNCOCIsCiJIIwljICNBMjk4QjciLAoiSSMJYyAjQkRCQUJBIiwKIkojCWMgI0I3QjVCNSIsCiJLIwljICNCMkIxQjEiLAoiTCMJYyAjQUNBREFEIiwKIk0jCWMgI0FCQURBQyIsCiJOIwljICNDQ0NEQ0MiLAoiTyMJYyAjRkJGQkZCIiwKIlAjCWMgIzg2NzNBRCIsCiJRIwljICM4MDcxQUEiLAoiUiMJYyAjN0E2RUE3IiwKIlMjCWMgIzgyNkZBOCIsCiJUIwljICM5MTcwQTgiLAoiVSMJYyAjNUU0NkE1IiwKIlYjCWMgI0U2QUJBQiIsCiJXIwljICNFNEFDQUIiLAoiWCMJYyAjREVBQ0FDIiwKIlkjCWMgI0Q3QUNBRCIsCiJaIwljICNENUFEQUUiLAoiYCMJYyAjNDUzN0E1IiwKIiAkCWMgIzM2MjhBNCIsCiIuJAljICM2RTUwQTYiLAoiKyQJYyAjOUI3MEE5IiwKIkAkCWMgI0JDODdBQSIsCiIjJAljICNERDlGQUMiLAoiJCQJYyAjRjFBRUFFIiwKIiUkCWMgI0Y0QUZBRiIsCiImJAljICNGOUIyQjMiLAoiKiQJYyAjRkNCNUI1IiwKIj0kCWMgI0ZGQkZCRiIsCiItJAljICNGRkMyQzIiLAoiOyQJYyAjRkZDNUM2IiwKIj4kCWMgIzg0NjhCNiIsCiIsJAljICNDREMwQzUiLAoiJyQJYyAjQzZDNEMzIiwKIikkCWMgI0MwQkVCRSIsCiIhJAljICNCOUI4QjgiLAoifiQJYyAjQjFCMUIyIiwKInskCWMgI0FEQUVBRCIsCiJdJAljICNBREFEQUQiLAoiXiQJYyAjRDJEMkQyIiwKIi8kCWMgI0UzQUNBRCIsCiIoJAljICNEREFFQUMiLAoiXyQJYyAjRDZCMkFEIiwKIjokCWMgI0NBQjBBQyIsCiI8JAljICNDNEIxQUMiLAoiWyQJYyAjQkZCMUFDIiwKIn0kCWMgI0JFQjBBRCIsCiJ8JAljICNCMEEyQUIiLAoiMSQJYyAjMkMyN0EzIiwKIjIkCWMgIzE4MTFBMiIsCiIzJAljICM1NTNEQTUiLAoiNCQJYyAjOTQ2QkFBIiwKIjUkCWMgI0UzQTJBRiIsCiI2JAljICNGQUIzQjQiLAoiNyQJYyAjRkZCOUJBIiwKIjgkCWMgI0ZGQkRCRCIsCiI5JAljICNGRkM2QzYiLAoiMCQJYyAjRkZDOUNBIiwKImEkCWMgIzZFNThCNSIsCiJiJAljICMwQTA4QTMiLAoiYyQJYyAjREJDRkNGIiwKImQkCWMgI0NGQ0RDRCIsCiJlJAljICNDOEM4QzciLAoiZiQJYyAjQzFDMUMxIiwKImckCWMgI0I4QjhCOCIsCiJoJAljICNCMkIyQjIiLAoiaSQJYyAjRUNFQ0VDIiwKImokCWMgI0U1QUJBQyIsCiJrJAljICNFMkFDQUYiLAoibCQJYyAjRDZBREFFIiwKIm0kCWMgI0M5QjJBRiIsCiJuJAljICNCRUIyQUQiLAoibyQJYyAjQjlCNkIwIiwKInAkCWMgI0IxQjJBRCIsCiJxJAljICNBREFGQUQiLAoiciQJYyAjQjBCMkFGIiwKInMkCWMgI0IwQjBCMCIsCiJ0JAljICNCMkIwQjIiLAoidSQJYyAjNzY3MUFCIiwKInYkCWMgIzE3MTVBMiIsCiJ3JAljICMxMTBEQTMiLAoieCQJYyAjOUM3MEFEIiwKInkkCWMgI0ZFQjdCOCIsCiJ6JAljICNGRkJBQkIiLAoiQSQJYyAjRkZDM0MzIiwKIkIkCWMgI0ZGQzhDOCIsCiJDJAljICNGRkNDQ0QiLAoiRCQJYyAjNjA0RkIzIiwKIkUkCWMgIzEzMTFBNiIsCiJGJAljICNFNEQ5RDkiLAoiRyQJYyAjRDlEOEQ3IiwKIkgkCWMgI0QyRDNEMSIsCiJJJAljICNDQ0NCQ0IiLAoiSiQJYyAjQzJDMUMyIiwKIkskCWMgI0I5QjlCOSIsCiJMJAljICNBQ0FDQUMiLAoiTSQJYyAjQkRCREJEIiwKIk4kCWMgI0Y0RjRGNCIsCiJPJAljICNERERERjIiLAoiUCQJYyAjMDMwM0ExIiwKIlEkCWMgI0I5QjlCMiIsCiJSJAljICNCN0I4QjQiLAoiUyQJYyAjQURBQ0FFIiwKIlQkCWMgI0I5QUJBQyIsCiJVJAljICNFMkFEQUQiLAoiViQJYyAjNUY0NkE1IiwKIlckCWMgI0U0QURBRiIsCiJYJAljICNENkFEQUYiLAoiWSQJYyAjQzVCMEFGIiwKIlokCWMgI0I4QjVCMiIsCiJgJAljICNCMkIyQjEiLAoiICUJYyAjQjJCM0IyIiwKIi4lCWMgI0IzQjRCMyIsCiIrJQljICNCNEI0QjQiLAoiQCUJYyAjQjVCNUI1IiwKIiMlCWMgI0JBQjlCQSIsCiIkJQljICNCQ0JCQkQiLAoiJSUJYyAjQUJBOUFDIiwKIiYlCWMgIzhEODVBOSIsCiIqJQljICM3MjVGQTgiLAoiPSUJYyAjNDgzNkE0IiwKIi0lCWMgIzFBMTNBMiIsCiI7JQljICM2MjQ3QUEiLAoiPiUJYyAjRkZCQkJDIiwKIiwlCWMgI0ZGQzBDMCIsCiInJQljICNGRkM1QzUiLAoiKSUJYyAjRkZDQUNBIiwKIiElCWMgI0ZGQ0VDRiIsCiJ+JQljICM3MDVFQjgiLAoieyUJYyAjMDQwNEEyIiwKIl0lCWMgI0VCRTFFMSIsCiJeJQljICNFMUUxRTAiLAoiLyUJYyAjREJEQ0RBIiwKIiglCWMgI0Q1RDVENSIsCiJfJQljICNDQ0NBQ0MiLAoiOiUJYyAjQzJDMEMwIiwKIjwlCWMgI0I2QjVCNiIsCiJbJQljICM4RDhEQUIiLAoifSUJYyAjNTQ1NEE3IiwKInwlCWMgIzYxNjFCQiIsCiIxJQljICM1QzVDQzMiLAoiMiUJYyAjQTlBOURGIiwKIjMlCWMgI0QyQ0VERiIsCiI0JQljICNDNEI3QjEiLAoiNSUJYyAjQkNCNkFGIiwKIjYlCWMgI0FFQUJBQiIsCiI3JQljICNCQ0FCQUMiLAoiOCUJYyAjRTVBREFFIiwKIjklCWMgIzM3MjhBMyIsCiIwJQljICNEOEFFQUYiLAoiYSUJYyAjQzJCNUIxIiwKImIlCWMgI0I0QjBBRSIsCiJjJQljICNBRkFGQUYiLAoiZCUJYyAjQjFCMkI0IiwKImUlCWMgI0IzQjNCNSIsCiJmJQljICNCMkIyQjQiLAoiZyUJYyAjQUZBRUFFIiwKImglCWMgI0FDQUJBQiIsCiJpJQljICNBQ0FDQUIiLAoiaiUJYyAjQkJBQkFCIiwKImslCWMgI0RDQUZCMiIsCiJsJQljICNFQUFDQUUiLAoibSUJYyAjRTdBNkFDIiwKIm4lCWMgI0M4OEZBOSIsCiJvJQljICM5NzZEQTkiLAoicCUJYyAjNTEzQUE1IiwKInElCWMgI0IwODJCNCIsCiJyJQljICNGRkMwQzEiLAoicyUJYyAjRkZDQkNDIiwKInQlCWMgI0ZGRDFEMiIsCiJ1JQljICM4NDcwQkUiLAoidiUJYyAjRTBEOEUzIiwKInclCWMgI0U5RTlFOSIsCiJ4JQljICNERkUwRTAiLAoieSUJYyAjRDVENUQ2IiwKInolCWMgI0NCQ0FDQSIsCiJBJQljICMzODM4QTciLAoiQiUJYyAjMjAyMEFDIiwKIkMlCWMgI0Y2RjNGNCIsCiJEJQljICNDNUJGQ0UiLAoiRSUJYyAjQzJCQUI1IiwKIkYlCWMgI0IzQUJBQiIsCiJHJQljICNDOEFEQUQiLAoiSCUJYyAjRTZBREFFIiwKIkklCWMgIzhENjdBNyIsCiJKJQljICMzNDI2QTMiLAoiSyUJYyAjQjg4N0E5IiwKIkwlCWMgI0RGQURBRiIsCiJNJQljICNDN0FEQUQiLAoiTiUJYyAjQjdCMUFEIiwKIk8lCWMgI0IzQjNCMiIsCiJQJQljICNBQ0FCQUMiLAoiUSUJYyAjQUJBQ0FDIiwKIlIlCWMgI0FCQUJBQiIsCiJTJQljICNCMUIxQjAiLAoiVCUJYyAjQURBQ0FCIiwKIlUlCWMgI0IzQUNBQiIsCiJWJQljICNDREFEQUYiLAoiVyUJYyAjRThBREIwIiwKIlglCWMgI0VFQURBQyIsCiJZJQljICNFRkFEQUIiLAoiWiUJYyAjRjFBREFEIiwKImAlCWMgI0YyQUZBRiIsCiIgJgljICNCOTg1QUMiLAoiLiYJYyAjMEEwN0EyIiwKIismCWMgIzRGM0FBOSIsCiJAJgljICNGRkNDQ0MiLAoiIyYJYyAjRkZEMkQyIiwKIiQmCWMgI0E4OEVDNSIsCiIlJgljICNBRkE5RDciLAoiJiYJYyAjRUZFRUVFIiwKIiomCWMgI0VDRUNFQSIsCiI9JgljICNFNkU3RTgiLAoiLSYJYyAjRERERURGIiwKIjsmCWMgI0Q0RDREMyIsCiI+JgljICNDNEM0QzQiLAoiLCYJYyAjNjY2NkFGIiwKIicmCWMgIzA2MDZBMyIsCiIpJgljICM1OTU5QzIiLAoiISYJYyAjRTdFM0UyIiwKIn4mCWMgI0I3QjBCRSIsCiJ7JgljICNCRkMxQjciLAoiXSYJYyAjQjBBRkFFIiwKIl4mCWMgI0JBQUJBQiIsCiIvJgljICNEOUFGQjAiLAoiKCYJYyAjRTlBQ0FEIiwKIl8mCWMgI0VCQUNBQiIsCiI6JgljICNFN0FCQUMiLAoiPCYJYyAjOUQ4MEFDIiwKIlsmCWMgIzA4MDhBMiIsCiJ9JgljICMyODI3QTUiLAoifCYJYyAjQjZCNkIzIiwKIjEmCWMgI0I3QjdCNCIsCiIyJgljICNCMkIwQjAiLAoiMyYJYyAjNzE3MUE4IiwKIjQmCWMgIzU2NTdBNyIsCiI1JgljICM5ODlBQUQiLAoiNiYJYyAjQjBCMkIxIiwKIjcmCWMgI0FFQUVBRCIsCiI4JgljICNCN0I2QUYiLAoiOSYJYyAjQzFBREFDIiwKIjAmCWMgI0VEQURBQiIsCiJhJgljICNGMUFGQUYiLAoiYiYJYyAjRjNCMEIwIiwKImMmCWMgIzdGNUNBQSIsCiJkJgljICMyMzFBQTUiLAoiZSYJYyAjRThDNUQ0IiwKImYmCWMgIzZFNkFDNCIsCiJnJgljICNGM0YyRjMiLAoiaCYJYyAjRjFGMkYxIiwKImkmCWMgI0VERUVFRSIsCiJqJgljICNFNEU2RTYiLAoiayYJYyAjREJEQ0RDIiwKImwmCWMgI0NFQ0VDRSIsCiJtJgljICNCQkJCQkIiLAoibiYJYyAjRjFGMUYxIiwKIm8mCWMgIzJCMkJCMSIsCiJwJgljICNCRkJGRTgiLAoicSYJYyAjRTBEQUQ5IiwKInImCWMgI0IxQThCNSIsCiJzJgljICNDOUM3QjkiLAoidCYJYyAjQzVBQ0FDIiwKInUmCWMgI0U0QURBRCIsCiJ2JgljICNFQ0FDQUMiLAoidyYJYyAjNjU1NUE4IiwKIngmCWMgI0FGQUZCOSIsCiJ5JgljICNCNUI2QjIiLAoieiYJYyAjQjJCMkIzIiwKIkEmCWMgI0I0QjJCMyIsCiJCJgljICNBNEE0QUEiLAoiQyYJYyAjMjYyN0E1IiwKIkQmCWMgI0IxQjJCMSIsCiJFJgljICNBRUFFQUMiLAoiRiYJYyAjQzBBREFDIiwKIkcmCWMgI0VDQURBQiIsCiJIJgljICNDMThDQUYiLAoiSSYJYyAjRjlCQkJGIiwKIkomCWMgI0ZGQ0JDQiIsCiJLJgljICNGRkQxRDEiLAoiTCYJYyAjRkZEOEQ4IiwKIk0mCWMgIzJFMjhBQyIsCiJOJgljICMxMTEwQTYiLAoiTyYJYyAjRUJFOUYyIiwKIlAmCWMgI0Y1RjVGNCIsCiJRJgljICNGMkYyRjEiLAoiUiYJYyAjRUFFQkVBIiwKIlMmCWMgI0UwRTJFMiIsCiJUJgljICNEM0Q0RDQiLAoiVSYJYyAjQzZDNkM2IiwKIlYmCWMgI0FFQUVBRSIsCiJXJgljICNEQkRCREIiLAoiWCYJYyAjRURFN0U2IiwKIlkmCWMgI0QzQzhDOCIsCiJaJgljICNBREE2QjQiLAoiYCYJYyAjQ0NCOEIyIiwKIiAqCWMgI0Q3QUVBRiIsCiIuKgljICNFREFDQUQiLAoiKyoJYyAjRTdBREFCIiwKIkAqCWMgI0UxQUVBRCIsCiIjKgljICM0RjQ1QTciLAoiJCoJYyAjOTg5OEFEIiwKIiUqCWMgI0I1QjVCMiIsCiImKgljICNCMEIxQjAiLAoiKioJYyAjQkRCREJBIiwKIj0qCWMgI0IxQjFBRSIsCiItKgljICM5MzkyQUMiLAoiOyoJYyAjOTk5N0IzIiwKIj4qCWMgI0FGQUVBRCIsCiIsKgljICNBRUFGQUQiLAoiJyoJYyAjQjhCOUIzIiwKIikqCWMgI0JFQURBRCIsCiIhKgljICNFMUFDQUYiLAoifioJYyAjRUNBQ0FEIiwKInsqCWMgI0YwQUNBRCIsCiJdKgljICNGMkFGQjAiLAoiXioJYyAjQkE4N0FFIiwKIi8qCWMgIzA3MDVBMiIsCiIoKgljICNGRkJFQkYiLAoiXyoJYyAjRkZDM0M0IiwKIjoqCWMgI0ZGRDBDRiIsCiI8KgljICNGRkQ3RDYiLAoiWyoJYyAjOUI4NkM2IiwKIn0qCWMgIzc4NzhDQSIsCiJ8KgljICNGNkY2RjYiLAoiMSoJYyAjRUZFRkVGIiwKIjIqCWMgI0U2RTdFNyIsCiIzKgljICNEOURBREEiLAoiNCoJYyAjQ0RDRENEIiwKIjUqCWMgI0JFQkVCRSIsCiI2KgljICNCM0IzQjMiLAoiNyoJYyAjRjlGOUY5IiwKIjgqCWMgI0VERTFFMSIsCiI5KgljICNDQ0JFQkYiLAoiMCoJYyAjQTVBOUI3IiwKImEqCWMgI0NBQjRCMCIsCiJiKgljICNDQkFDQUIiLAoiYyoJYyAjRTRBRUFFIiwKImQqCWMgI0VFQUJBQiIsCiJlKgljICNFMkFFQUUiLAoiZioJYyAjNEU0NUE3IiwKImcqCWMgIzkyOTJBQSIsCiJoKgljICNCM0IzQjEiLAoiaSoJYyAjQkRCREI3IiwKImoqCWMgI0JFQkRCNiIsCiJrKgljICNBMkExQjEiLAoibCoJYyAjNEY0REFCIiwKIm0qCWMgI0IwQUVBRCIsCiJuKgljICNCNUI3QjEiLAoibyoJYyAjQkZBREFEIiwKInAqCWMgI0U1QURCMCIsCiJxKgljICM3RDVBQTkiLAoicioJYyAjMzIyNEE1IiwKInMqCWMgI0ZGQkRCRSIsCiJ0KgljICNGRkMyQzMiLAoidSoJYyAjRkZDN0M4IiwKInYqCWMgI0ZGQ0VDRCIsCiJ3KgljICNGRkQ1RDQiLAoieCoJYyAjRjhENkRBIiwKInkqCWMgIzFDMTlBOCIsCiJ6KgljICNCNkI2RTAiLAoiQSoJYyAjRjVGNUY1IiwKIkIqCWMgI0YxRjBGMCIsCiJDKgljICNFOUVBRUEiLAoiRCoJYyAjREREREREIiwKIkUqCWMgI0QxRDFEMSIsCiJGKgljICNCNkI2QjYiLAoiRyoJYyAjRThFOEU4IiwKIkgqCWMgI0M4QzhDOCIsCiJJKgljICNDN0M3QzciLAoiSioJYyAjRkFGOEY4IiwKIksqCWMgI0U3RDNENCIsCiJMKgljICNDQ0I3QjgiLAoiTSoJYyAjQTdCM0MwIiwKIk4qCWMgI0M0QUJBQyIsCiJPKgljICNFMkFFQUQiLAoiUCoJYyAjRUZBQkFCIiwKIlEqCWMgI0U1QURBRCIsCiJSKgljICM0RjQ1QTYiLAoiUyoJYyAjOTE5MkFBIiwKIlQqCWMgI0FEQURBQyIsCiJVKgljICNCMkIxQUYiLAoiVioJYyAjQkNCQkI2IiwKIlcqCWMgI0FCQUJCMiIsCiJYKgljICM3QTc5QTgiLAoiWSoJYyAjQkNCN0IxIiwKIloqCWMgI0NBQUVBRCIsCiJgKgljICNFOUFFQUYiLAoiID0JYyAjQzM4REFDIiwKIi49CWMgIzg4NjNBRSIsCiIrPQljICNGRUJDQkMiLAoiQD0JYyAjRkZEMkQxIiwKIiM9CWMgI0ZFRDhENyIsCiIkPQljICNCREE2Q0YiLAoiJT0JYyAjMDkwOUE0IiwKIiY9CWMgI0E3QTZEOSIsCiIqPQljICNGMEVGRjAiLAoiPT0JYyAjRUFFQUVBIiwKIi09CWMgI0RGREZERiIsCiI7PQljICNENEQ0RDQiLAoiPj0JYyAjQzVDNUM1IiwKIiw9CWMgI0RFREVERSIsCiInPQljICNENkQ2RDYiLAoiKT0JYyAjRDhEOEQ4IiwKIiE9CWMgI0NBQ0FDQSIsCiJ+PQljICNCNUIzREUiLAoiez0JYyAjNTk0REIwIiwKIl09CWMgIzNBMzFBNSIsCiJePQljICMyRDJFQTkiLAoiLz0JYyAjMzYyQkE1IiwKIig9CWMgIzNCMkJBNSIsCiJfPQljICMzQzJCQTQiLAoiOj0JYyAjM0IyQkE0IiwKIjw9CWMgIzNBMkJBNCIsCiJbPQljICMxNDExQTIiLAoifT0JYyAjOTE5MUFCIiwKInw9CWMgI0FEQURBRiIsCiIxPQljICNBRUFFQjAiLAoiMj0JYyAjQUJBQkFDIiwKIjM9CWMgI0IwQjBCMSIsCiI0PQljICMwMjAyQTEiLAoiNT0JYyAjNTI1MUE2IiwKIjY9CWMgI0FDOURBQyIsCiI3PQljICNEOUIwQUMiLAoiOD0JYyAjRTlBQ0FDIiwKIjk9CWMgI0U2QTdBQyIsCiIwPQljICM4NDYwQTgiLAoiYT0JYyAjMTIwREEyIiwKImI9CWMgIzFFMTZBMyIsCiJjPQljICNGMkIwQjYiLAoiZD0JYyAjRkRCQkJCIiwKImU9CWMgI0ZFQkZCRiIsCiJmPQljICNGRkM4QzkiLAoiZz0JYyAjRkZDRkNGIiwKImg9CWMgI0ZGRDVENSIsCiJpPQljICNGRkRDREMiLAoiaj0JYyAjOEE3QUM0IiwKIms9CWMgIzRDNENCQSIsCiJsPQljICNCMkIyRDgiLAoibT0JYyAjRERERERFIiwKIm49CWMgI0E4QThBRiIsCiJvPQljICM2RTZFQTgiLAoicD0JYyAjMkUyRUE5IiwKInE9CWMgI0RBREFEQSIsCiJyPQljICM1NTUzQkIiLAoicz0JYyAjOTY5N0FDIiwKInQ9CWMgI0FFQURBRSIsCiJ1PQljICMwNzA3QTIiLAoidj0JYyAjMkIyMUEzIiwKInc9CWMgIzQ3MzRBNCIsCiJ4PQljICM0RDM4QTUiLAoieT0JYyAjMkUyMUEzIiwKIno9CWMgIzBBMDdBMSIsCiJBPQljICMxMzBFQTIiLAoiQj0JYyAjQ0M5NEFGIiwKIkM9CWMgI0ZBQjVCNSIsCiJEPQljICNGQ0I5QjkiLAoiRT0JYyAjRkVCREJEIiwKIkY9CWMgI0ZGRDdENyIsCiJHPQljICNGRkRDREQiLAoiSD0JYyAjOTQ4OUM5IiwKIkk9CWMgIzEwMTBBNiIsCiJKPQljICMyQTJBQUIiLAoiSz0JYyAjMzAzMEFBIiwKIkw9CWMgIzE2MTZBNCIsCiJNPQljICNDNUM1RTMiLAoiTj0JYyAjRDNEM0QzIiwKIk89CWMgIzlDOTZDOSIsCiJQPQljICMxNDEwQTMiLAoiUT0JYyAjMzgzNkE0IiwKIlI9CWMgI0FGQUZBRSIsCiJTPQljICNCMEFGQUQiLAoiVD0JYyAjQjNCMkIwIiwKIlU9CWMgIzFGMUZBNCIsCiJWPQljICMyMjIyQTUiLAoiVz0JYyAjNDQzMUE1IiwKIlg9CWMgI0Q2OUFBRSIsCiJZPQljICNGNUIxQjEiLAoiWj0JYyAjRjdCNEI0IiwKImA9CWMgI0ZBQjdCNyIsCiIgLQljICNGQ0JBQkEiLAoiLi0JYyAjRkVCRUJFIiwKIistCWMgI0ZGQzdDOSIsCiJALQljICNGRkNEQ0QiLAoiIy0JYyAjRkVEMkQzIiwKIiQtCWMgI0ZFRDhEOCIsCiIlLQljICNGMkREREQiLAoiJi0JYyAjQkNCOUQ1IiwKIiotCWMgIzQyNDJCNSIsCiI9LQljICMyNjI2QUYiLAoiLS0JYyAjN0I3QkMwIiwKIjstCWMgI0MwQzBDMCIsCiI+LQljICNFNEQzRDMiLAoiLC0JYyAjRUFCOEI4IiwKIictCWMgI0Y2QUNBQyIsCiIpLQljICNGNUFEQUUiLAoiIS0JYyAjRjRBQ0FEIiwKIn4tCWMgI0Y1QUJBRCIsCiJ7LQljICNGNEFCQUMiLAoiXS0JYyAjRjJBQkFCIiwKIl4tCWMgI0YxQUJBRCIsCiIvLQljICNFQUFDQUMiLAoiKC0JYyAjQzdBQ0FDIiwKIl8tCWMgI0I1QjJCMSIsCiI6LQljICNCMkIzQjEiLAoiPC0JYyAjQUZBREFFIiwKIlstCWMgIzU3NTdBNyIsCiJ9LQljICMwQjBCQTIiLAoifC0JYyAjQjRBOUIyIiwKIjEtCWMgI0I0OUZBQyIsCiIyLQljICM4MDY1QTkiLAoiMy0JYyAjNDYzNUE1IiwKIjQtCWMgIzJCMjBBMyIsCiI1LQljICMwRTBCQTIiLAoiNi0JYyAjMDYwNUEyIiwKIjctCWMgIzE1MTBBMiIsCiI4LQljICM2RDUwQTYiLAoiOS0JYyAjQzE4Q0FDIiwKIjAtCWMgI0VEQUVBRiIsCiJhLQljICNFRkFGQUYiLAoiYi0JYyAjRjFCMEIxIiwKImMtCWMgI0YzQjNCMiIsCiJkLQljICNGNkI1QjUiLAoiZS0JYyAjRjlCOEI4IiwKImYtCWMgI0ZCQkJCQyIsCiJnLQljICNGQkJGQzAiLAoiaC0JYyAjRkJDM0M0IiwKImktCWMgI0ZCQzhDOCIsCiJqLQljICNGQkNEQ0QiLAoiay0JYyAjRkJEMkQzIiwKImwtCWMgI0VERDZENiIsCiJtLQljICNEQ0Q5RDgiLAoibi0JYyAjRENEQ0RDIiwKIm8tCWMgI0I1QjVEMiIsCiJwLQljICM3MjcyQzAiLAoicS0JYyAjMzYzN0FGIiwKInItCWMgIzFDMUNBOCIsCiJzLQljICMwOTA5QTMiLAoidC0JYyAjMDUwNUEyIiwKInUtCWMgIzEwMTBBMyIsCiJ2LQljICMxRDFEQTQiLAoidy0JYyAjMzkzOUE1IiwKIngtCWMgIzVCNUJBOCIsCiJ5LQljICNBREFEQzciLAoiei0JYyAjRTRFNEU0IiwKIkEtCWMgI0U5RDBEMSIsCiJCLQljICNFQUIzQjQiLAoiQy0JYyAjRjdBQkFDIiwKIkQtCWMgI0ZBQUJBQiIsCiJFLQljICNGOUFCQUMiLAoiRi0JYyAjRjhBQkFEIiwKIkctCWMgI0Y2QUJBQyIsCiJILQljICNGNUFCQUIiLAoiSS0JYyAjRjRBQkFCIiwKIkotCWMgI0YwQUJBQiIsCiJLLQljICNFMUFEQjAiLAoiTC0JYyAjQ0ZBRUFGIiwKIk0tCWMgI0I5QUJBQiIsCiJOLQljICNCMUFDQUIiLAoiTy0JYyAjQjFBQ0FEIiwKIlAtCWMgI0IyQUNBRCIsCiJRLQljICNCMUFEQUQiLAoiUi0JYyAjOUE5NEFCIiwKIlMtCWMgI0I2QThBQyIsCiJULQljICNDRkIwQUYiLAoiVS0JYyAjREJBREFDIiwKIlYtCWMgI0U2QUJBRCIsCiJXLQljICNFN0FDQUQiLAoiWC0JYyAjRTlBREFFIiwKIlktCWMgI0VBQUVBRiIsCiJaLQljICNFQ0IwQjAiLAoiYC0JYyAjRUVCMkIxIiwKIiA7CWMgI0YxQjRCMyIsCiIuOwljICNGNEI2QjYiLAoiKzsJYyAjRjZCOUJBIiwKIkA7CWMgI0Y2QkNCRCIsCiIjOwljICNGNkMwQzAiLAoiJDsJYyAjRjdDNEM0IiwKIiU7CWMgI0Y3QzhDOCIsCiImOwljICNGN0NDQ0QiLAoiKjsJYyAjRTdEMEQwIiwKIj07CWMgI0Q0RDFEMCIsCiItOwljICNENEQzRDQiLAoiOzsJYyAjQ0ZDRkNGIiwKIj47CWMgI0MyQzJDMiIsCiIsOwljICNCQUJBQkEiLAoiJzsJYyAjRURDRUNFIiwKIik7CWMgI0VEQjFCMyIsCiIhOwljICNGOEFDQUUiLAoifjsJYyAjRjlBQkFCIiwKIns7CWMgI0Y4QUJBQyIsCiJdOwljICNGN0FDQUIiLAoiXjsJYyAjRjZBQ0FCIiwKIi87CWMgI0YzQUJBQyIsCiIoOwljICNGMEFCQUQiLAoiXzsJYyAjRTdBQ0FDIiwKIjo7CWMgI0REQUZBRSIsCiI8OwljICNEM0FDQUMiLAoiWzsJYyAjQzlBQ0FEIiwKIn07CWMgI0M3QURBQyIsCiJ8OwljICNDQ0FDQUMiLAoiMTsJYyAjRDVBREFEIiwKIjI7CWMgI0RBQUZBRSIsCiIzOwljICNEREFEQUMiLAoiNDsJYyAjRERBQ0FDIiwKIjU7CWMgI0RDQUJBQyIsCiI2OwljICNEREFCQUIiLAoiNzsJYyAjREZBQkFCIiwKIjg7CWMgI0UzQUVBRSIsCiI5OwljICNFNUFGQUYiLAoiMDsJYyAjRTdCMEIwIiwKImE7CWMgI0U5QjJCMSIsCiJiOwljICNFQkIzQjMiLAoiYzsJYyAjRUVCNkI2IiwKImQ7CWMgI0VGQjlCOSIsCiJlOwljICNGMEJDQkMiLAoiZjsJYyAjRjBCRkJGIiwKImc7CWMgI0YwQzJDMiIsCiJoOwljICNGMEM1QzYiLAoiaTsJYyAjREZDOEM4IiwKImo7CWMgI0NBQzdDNyIsCiJrOwljICNDOUM4QzkiLAoibDsJYyAjQjdCN0I3IiwKIm07CWMgI0Q5RDlEOSIsCiJuOwljICNFRENDQ0QiLAoibzsJYyAjRjhBQ0FGIiwKInA7CWMgI0Y4QUNBRCIsCiJxOwljICNGN0FCQUIiLAoicjsJYyAjRjNBQ0FCIiwKInM7CWMgI0VCQUNBRCIsCiJ0OwljICNFQkFDQUMiLAoidTsJYyAjRThBQ0FDIiwKInY7CWMgI0UyQURBRSIsCiJ3OwljICNERkFFQjAiLAoieDsJYyAjRERBRkFGIiwKInk7CWMgI0RGQURBRCIsCiJ6OwljICNEREFDQUQiLAoiQTsJYyAjREFBQ0FDIiwKIkI7CWMgI0Q4QUJBQyIsCiJDOwljICNENkFCQUMiLAoiRDsJYyAjRDdBQ0FCIiwKIkU7CWMgI0Q4QUNBQiIsCiJGOwljICNEOUFDQUMiLAoiRzsJYyAjRDlBQ0FEIiwKIkg7CWMgI0Q5QURBRCIsCiJJOwljICNEQkFEQUUiLAoiSjsJYyAjRENBRUFFIiwKIks7CWMgI0RGQUZCMCIsCiJMOwljICNFMUIxQjEiLAoiTTsJYyAjRTJCM0IyIiwKIk47CWMgI0U1QjVCNSIsCiJPOwljICNFNkI3QjciLAoiUDsJYyAjRTdCOUJBIiwKIlE7CWMgI0U3QkNCQyIsCiJSOwljICNFN0JGQkYiLAoiUzsJYyAjRTdDMEMwIiwKIlQ7CWMgI0Q4QzJDMiIsCiJVOwljICNDM0JGQzAiLAoiVjsJYyAjQzFDMEMxIiwKIlc7CWMgI0MxQzBDMCIsCiJYOwljICNCRkJGQkYiLAoiWTsJYyAjQjFCMUIxIiwKIlo7CWMgI0VCQ0NDQyIsCiJgOwljICNFOEFGQjAiLAoiID4JYyAjREJBQkFCIiwKIi4+CWMgI0Q5QUJBQiIsCiIrPgljICNENUFCQUMiLAoiQD4JYyAjRDJBQkFDIiwKIiM+CWMgI0NGQUJBQiIsCiIkPgljICNDREFCQUIiLAoiJT4JYyAjQ0NBQkFCIiwKIiY+CWMgI0NFQUJBQiIsCiIqPgljICNEMEFCQUIiLAoiPT4JYyAjRDJBQkFCIiwKIi0+CWMgI0QzQUJBQiIsCiI7PgljICNEMkFDQUIiLAoiPj4JYyAjRDVBQkFEIiwKIiw+CWMgI0Q1QUNBRCIsCiInPgljICNENkFFQUQiLAoiKT4JYyAjRDhBRkIwIiwKIiE+CWMgI0Q5QjBCMSIsCiJ+PgljICNEQUIxQjIiLAoiez4JYyAjREFCM0I0IiwKIl0+CWMgI0RBQjVCNCIsCiJePgljICNEQUI2QjUiLAoiLz4JYyAjQ0VCOUI5IiwKIig+CWMgI0JFQjlCQSIsCiJfPgljICNCQkI5QkEiLAoiOj4JYyAjQkJCOUI5IiwKIjw+CWMgI0JBQkJCQiIsCiJbPgljICNFQUNEQ0MiLAoifT4JYyAjRTZCMkIxIiwKInw+CWMgI0VEQUVBRCIsCiIxPgljICNFQkFFQUUiLAoiMj4JYyAjRTlBRUFFIiwKIjM+CWMgI0VBQUVBRSIsCiI0PgljICNFQUFEQUQiLAoiNT4JYyAjRThBRUFFIiwKIjY+CWMgI0U3QUVBRCIsCiI3PgljICNFNkFFQUQiLAoiOD4JYyAjRTZBREFGIiwKIjk+CWMgI0U1QURBRiIsCiIwPgljICNFM0FEQUUiLAoiYT4JYyAjRTNBREFDIiwKImI+CWMgI0UxQURBQyIsCiJjPgljICNERkFEQUMiLAoiZD4JYyAjREVBREFDIiwKImU+CWMgI0RDQURBQyIsCiJmPgljICNEQkFEQUQiLAoiZz4JYyAjREJBQ0FDIiwKImg+CWMgI0Q1QUNBQyIsCiJpPgljICNEMUFDQUMiLAoiaj4JYyAjQ0NBQ0FCIiwKIms+CWMgI0NBQUNBQiIsCiJsPgljICNDQUFDQUMiLAoibT4JYyAjQ0RBQ0FDIiwKIm4+CWMgI0NFQUNBQyIsCiJvPgljICNDRkFDQUIiLAoicD4JYyAjQ0ZBQ0FDIiwKInE+CWMgI0QxQUJBQyIsCiJyPgljICNEMkFDQUMiLAoicz4JYyAjRDNBQ0FEIiwKInQ+CWMgI0Q0QUNBQyIsCiJ1PgljICNENEFEQUIiLAoidj4JYyAjRDVBRUFGIiwKInc+CWMgI0RBQjNCMiIsCiJ4PgljICNDREI2QjYiLAoieT4JYyAjQkFCNkI3IiwKIno+CWMgI0I4QjdCNyIsCiJBPgljICNCOEI3QjYiLAoiQj4JYyAjQjdCN0I4IiwKIkM+CWMgI0VEQ0RDRCIsCiJEPgljICNFRUIyQjIiLAoiRT4JYyAjRkFBREFEIiwKIkY+CWMgI0Y5QURBQyIsCiJHPgljICNGN0FDQUUiLAoiSD4JYyAjRjZBREFDIiwKIkk+CWMgI0Y1QUNBQiIsCiJKPgljICNGNUFDQUMiLAoiSz4JYyAjRjFBQ0FCIiwKIkw+CWMgI0YxQUNBRCIsCiJNPgljICNGMEFCQUUiLAoiTj4JYyAjRTZBQ0FEIiwKIk8+CWMgI0U1QUNBRCIsCiJQPgljICNENkFDQUIiLAoiUT4JYyAjRDBBQ0FDIiwKIlI+CWMgI0QwQUJBQyIsCiJTPgljICNENEFCQUIiLAoiVD4JYyAjRDRBQkFDIiwKIlU+CWMgI0Q1QUJBQiIsCiJWPgljICNENkFCQUIiLAoiVz4JYyAjRDdBQkFEIiwKIlg+CWMgI0Q4QURBQiIsCiJZPgljICNEOEFFQUQiLAoiWj4JYyAjRENCMEIwIiwKImA+CWMgI0REQjBCMCIsCiIgLAljICNERUIxQjEiLAoiLiwJYyAjREVCMkIyIiwKIissCWMgI0NGQjVCNSIsCiJALAljICNCOUI0QjQiLAoiIywJYyAjQjVCNEI0IiwKIiQsCWMgI0I1QjVCNCIsCiIlLAljICNEN0Q3RDciLAoiJiwJYyAjRjBCMUIxIiwKIiosCWMgI0ZDQUJBQiIsCiI9LAljICNGQ0FCQUMiLAoiLSwJYyAjRkJBQkFDIiwKIjssCWMgI0Y1QUJBQyIsCiI+LAljICNGMUFCQUIiLAoiLCwJYyAjRERBQ0FCIiwKIicsCWMgI0Q3QUJBQiIsCiIpLAljICNEOEFCQUIiLAoiISwJYyAjREFBQkFDIiwKIn4sCWMgI0Q5QUNBQiIsCiJ7LAljICNEQUFEQUIiLAoiXSwJYyAjRENBREFCIiwKIl4sCWMgI0REQURBRCIsCiIvLAljICNERkFFQUYiLAoiKCwJYyAjRTBBRkFGIiwKIl8sCWMgI0UyQjBCMCIsCiI6LAljICNFM0IxQjEiLAoiPCwJYyAjRDFCNEI0IiwKIlssCWMgI0I4QjJCMyIsCiJ9LAljICNCNEIzQjMiLAoifCwJYyAjQjRCNEIzIiwKIjEsCWMgI0I0QjNCNCIsCiIyLAljICNFQUQxRDAiLAoiMywJYyAjRUVCNUIzIiwKIjQsCWMgI0ZDQUNBQyIsCiI1LAljICNGQkFCQUIiLAoiNiwJYyAjRjhBQkFCIiwKIjcsCWMgI0YzQUJBQiIsCiI4LAljICNEQkFCQUMiLAoiOSwJYyAjRERBQkFEIiwKIjAsCWMgI0RFQUJBRCIsCiJhLAljICNERUFEQUIiLAoiYiwJYyAjRTBBREFCIiwKImMsCWMgI0UyQUNBRCIsCiJkLAljICNFM0FDQUUiLAoiZSwJYyAjRTRBRUFGIiwKImYsCWMgI0U2QUVBRiIsCiJnLAljICNFOEIxQjEiLAoiaCwJYyAjRDRCNEI0IiwKImksCWMgI0I4QjFCMiIsCiJqLAljICNCM0IyQjMiLAoiaywJYyAjQjNCMkIyIiwKImwsCWMgI0U2RDREMiIsCiJtLAljICNFQ0I5QjciLAoibiwJYyAjRkFBQ0FDIiwKIm8sCWMgI0ZDQUJBRiIsCiJwLAljICNGOEFEQjAiLAoicSwJYyAjRjVBREFGIiwKInIsCWMgI0Y0QURBRCIsCiJzLAljICNFMkFCQUUiLAoidCwJYyAjRTFBQkFEIiwKInUsCWMgI0RFQURBRCIsCiJ2LAljICNERkFEQUIiLAoidywJYyAjRTFBREFCIiwKIngsCWMgI0U4QUNBRiIsCiJ5LAljICNFN0FEQUYiLAoieiwJYyAjRThBREFFIiwKIkEsCWMgI0VDQjBCMSIsCiJCLAljICNENUIzQjMiLAoiQywJYyAjQjdCMEIwIiwKIkQsCWMgI0YyRjJGMiIsCiJFLAljICNFMEUwRTAiLAoiRiwJYyAjRjVFREVCIiwKIkcsCWMgI0VDQkRCQiIsCiJILAljICNGQUFDQUYiLAoiSSwJYyAjRjRBRkI1IiwKIkosCWMgI0Q5QUZCNCIsCiJLLAljICNENEI3QjkiLAoiTCwJYyAjRTZCM0IyIiwKIk0sCWMgI0Y2QUJBQiIsCiJOLAljICNGM0FCQUQiLAoiTywJYyAjRTNBQkFEIiwKIlAsCWMgI0UzQUNBQiIsCiJRLAljICNFMUFEQUQiLAoiUiwJYyAjREFBRUFFIiwKIlMsCWMgI0Q0QUZBRiIsCiJULAljICNEMUFGQUUiLAoiVSwJYyAjRDFBRkFEIiwKIlYsCWMgI0QyQUZBRCIsCiJXLAljICNEOUFFQjEiLAoiWCwJYyAjRDdBNEFFIiwKIlksCWMgIzMyMjVBMyIsCiJaLAljICNBRjgxQUEiLAoiYCwJYyAjRTlBQ0FFIiwKIiAnCWMgI0VBQURBRSIsCiIuJwljICNFQ0FEQUUiLAoiKycJYyAjRUVBREFFIiwKIkAnCWMgI0VFQUVBRiIsCiIjJwljICNENUIxQjEiLAoiJCcJYyAjQjVBRkFFIiwKIiUnCWMgI0IwQjBBRiIsCiImJwljICNCQ0JDQkMiLAoiKicJYyAjRjhGM0YzIiwKIj0nCWMgI0VBQkZCRiIsCiItJwljICNGNkFFQjIiLAoiOycJYyAjRTlCMkI5IiwKIj4nCWMgI0NCQjFCQiIsCiIsJwljICNCREIxQjkiLAoiJycJYyAjQzRCM0I2IiwKIiknCWMgI0Q0QjVCNSIsCiIhJwljICNEREFGQjAiLAoificJYyAjRUVBRkFFIiwKInsnCWMgI0Q4QjJBRSIsCiJdJwljICNDQkIzQjMiLAoiXicJYyAjQkRBQ0FFIiwKIi8nCWMgI0I4QUJBQiIsCiIoJwljICNCNkFCQUMiLAoiXycJYyAjQjZBQ0FEIiwKIjonCWMgI0MzQjZCOSIsCiI8JwljICNBOTkzQjUiLAoiWycJYyAjNjQ0QkE3IiwKIn0nCWMgI0VBQUNBRCIsCiJ8JwljICNFRUFEQUQiLAoiMScJYyAjRURBRUFFIiwKIjInCWMgI0Q0QjBCMSIsCiIzJwljICNCMkFEQUQiLAoiNCcJYyAjQUVBRUFGIiwKIjUnCWMgI0FFQUZBRiIsCiI2JwljICNGQUY0RjQiLAoiNycJYyAjRThDNEM0IiwKIjgnCWMgI0VFQjJCNSIsCiI5JwljICNERUIzQkIiLAoiMCcJYyAjQ0NCNkNBIiwKImEnCWMgI0NEQkJDRCIsCiJiJwljICNDRkI3QkMiLAoiYycJYyAjRDRCOUI1IiwKImQnCWMgI0M4QUZBRCIsCiJlJwljICNERUIyQUUiLAoiZicJYyAjRjJBREFCIiwKImcnCWMgI0U2QUNBQiIsCiJoJwljICM3QTVDQTciLAoiaScJYyAjQjU5NkIzIiwKImonCWMgI0RFQzdDMCIsCiJrJwljICNEMUNGQ0MiLAoibCcJYyAjQjlCQ0JEIiwKIm0nCWMgI0FDQUJBRSIsCiJuJwljICNBREFCQUMiLAoibycJYyAjQUNBQkFEIiwKInAnCWMgI0FGQjBCNSIsCiJxJwljICNDNkNFRDEiLAoicicJYyAjOEQ4N0I1IiwKInMnCWMgIzU3NDJBNSIsCiJ0JwljICNFQkFEQUQiLAoidScJYyAjRDJCMEIwIiwKInYnCWMgI0IwQUNBQyIsCiJ3JwljICNBREFEQUUiLAoieCcJYyAjQURBQ0FEIiwKInknCWMgI0YzRjNGMyIsCiJ6JwljICNGOEY0RjUiLAoiQScJYyAjRThDRENEIiwKIkInCWMgI0U2QjhCOSIsCiJDJwljICNENkIyQkIiLAoiRCcJYyAjQzNCM0NFIiwKIkUnCWMgI0NCQjdDRiIsCiJGJwljICNERUJBQkQiLAoiRycJYyAjRUFDM0I2IiwKIkgnCWMgI0UxQzJCNiIsCiJJJwljICNGMUNCQkMiLAoiSicJYyAjRjFCQUFFIiwKIksnCWMgI0YyQUNBQiIsCiJMJwljICM1RjQ3QTUiLAoiTScJYyAjOUY5M0JEIiwKIk4nCWMgI0QwQkRCMyIsCiJPJwljICNDRUNBQzQiLAoiUCcJYyAjQjhCQUJBIiwKIlEnCWMgI0FCQUNBRSIsCiJSJwljICNBQkFCQUQiLAoiUycJYyAjQUNBRUIxIiwKIlQnCWMgI0JFQzVDNyIsCiJVJwljICNCQkIyQkMiLAoiVicJYyAjMUUxREE2IiwKIlcnCWMgIzA2MDZBMiIsCiJYJwljICM4NTZCQTkiLAoiWScJYyAjRDBBRkFGIiwKIlonCWMgI0FGQUNBQiIsCiJgJwljICMzMDMwQTQiLAoiICkJYyAjMTQxNEEyIiwKIi4pCWMgI0E0QTRBQiIsCiIrKQljICNGMERGREYiLAoiQCkJYyAjRTJCRUJFIiwKIiMpCWMgI0Q4QjNCOCIsCiIkKQljICNDOEI2Q0UiLAoiJSkJYyAjQzlCNUNCIiwKIiYpCWMgI0RBQjVCQiIsCiIqKQljICNFN0JCQjQiLAoiPSkJYyAjRTZDMUI1IiwKIi0pCWMgI0Y1RDVCRSIsCiI7KQljICNGQUQxQjgiLAoiPikJYyAjRjFCMUFDIiwKIiwpCWMgI0RFQTZBQSIsCiInKQljICNDQTk5QUIiLAoiKSkJYyAjMDgwN0EyIiwKIiEpCWMgIzkxOTlDMCIsCiJ+KQljICNDQ0Q4QzgiLAoieykJYyAjQjlDM0M4IiwKIl0pCWMgI0I1QjlDMCIsCiJeKQljICNCREMxQzIiLAoiLykJYyAjQjFCN0JEIiwKIigpCWMgI0FEQjFCOCIsCiJfKQljICNCRkMyQzIiLAoiOikJYyAjQzZDOUM3IiwKIjwpCWMgI0I2QkZDNSIsCiJbKQljICNDOUQwQzgiLAoifSkJYyAjQ0RENUQwIiwKInwpCWMgI0IzQjVCOSIsCiIxKQljICNDNEFDQUMiLAoiMikJYyAjRTdBQ0FCIiwKIjMpCWMgI0NFQUVBRiIsCiI0KQljICNBRkFCQUIiLAoiNSkJYyAjQTNBMkFCIiwKIjYpCWMgIzhEOERBOSIsCiI3KQljICNGOEY4RjgiLAoiOCkJYyAjRjRFRkVGIiwKIjkpCWMgI0U0QzhDOCIsCiIwKQljICNEREIwQjUiLAoiYSkJYyAjQ0NCNUNBIiwKImIpCWMgI0NDQjhDQyIsCiJjKQljICNEQUI2QkMiLAoiZCkJYyAjRThCQUI1IiwKImUpCWMgI0U3QkZCMSIsCiJmKQljICNGM0QzQjYiLAoiZykJYyAjRjlENUIyIiwKImgpCWMgI0Y2QzBCMCIsCiJpKQljICNGMkFFQUIiLAoiaikJYyAjMDkwOEExIiwKImspCWMgIzkzOTdCRCIsCiJsKQljICNDREQ1QzQiLAoibSkJYyAjQkNDQ0NGIiwKIm4pCWMgI0M4Q0NDRCIsCiJvKQljICNDOUM1QzEiLAoicCkJYyAjQzJDQkNGIiwKInEpCWMgI0I0QkVDOCIsCiJyKQljICNDQ0NEQ0QiLAoicykJYyAjQzRDNEJGIiwKInQpCWMgI0M5RDZEQSIsCiJ1KQljICNDQkQ1Q0QiLAoidikJYyAjQzJDQUMxIiwKIncpCWMgI0JGQzRDNCIsCiJ4KQljICNEN0FFQUUiLAoieSkJYyAjQ0RBRUFFIiwKInopCWMgI0ExQTFBQSIsCiJBKQljICM4MTgxQTkiLAoiQikJYyAjRThENUQ0IiwKIkMpCWMgI0U1QjZCQyIsCiJEKQljICNDQ0IxQzMiLAoiRSkJYyAjQ0NCN0NBIiwKIkYpCWMgI0RDQjhDMCIsCiJHKQljICNFN0JBQjIiLAoiSCkJYyAjRTdCRUFGIiwKIkkpCWMgI0YzRDFCMiIsCiJKKQljICNGOUQ3QUUiLAoiSykJYyAjRjhDREFGIiwKIkwpCWMgI0YwQjNBQiIsCiJNKQljICNEQUFEQUUiLAoiTikJYyAjOTI4RkI2IiwKIk8pCWMgI0NFQzZCRCIsCiJQKQljICNDMUNEQ0IiLAoiUSkJYyAjQ0JDRUNBIiwKIlIpCWMgI0QxQzhDMiIsCiJTKQljICNCREMwQzEiLAoiVCkJYyAjQUZCM0I5IiwKIlUpCWMgI0M3QzZDNyIsCiJWKQljICNDQkNBQzciLAoiVykJYyAjQjhDMUMyIiwKIlgpCWMgI0M2Q0FDNCIsCiJZKQljICNDQ0NEQzciLAoiWikJYyAjQkZDMkMwIiwKImApCWMgI0IxQUJBQiIsCiIgIQljICNEMkFEQUQiLAoiLiEJYyAjQ0FBREFFIiwKIishCWMgI0EyQTFBQSIsCiJAIQljICNGNkVDRUMiLAoiIyEJYyAjRTVDMUM1IiwKIiQhCWMgI0NFQjJCRCIsCiIlIQljICNBMThFQkIiLAoiJiEJYyAjNTk0Q0FEIiwKIiohCWMgIzI3MUZBMyIsCiI9IQljICMwRTBCQTMiLAoiLSEJYyAjMDcwNkEyIiwKIjshCWMgIzIwMUNBMyIsCiI+IQljICM1QTREQTUiLAoiLCEJYyAjQjA4OUE5IiwKIichCWMgI0U0QTdBQiIsCiIpIQljICM4ODY0QTciLAoiISEJYyAjMzkyQUEzIiwKIn4hCWMgIzEyMERBMSIsCiJ7IQljICMxRDE1QTIiLAoiXSEJYyAjNTc0MUE1IiwKIl4hCWMgI0I4OEVBQSIsCiIvIQljICNEMkFDQUQiLAoiKCEJYyAjOTE5NEMxIiwKIl8hCWMgI0MzQzVDNiIsCiI6IQljICNBRUIyQjkiLAoiPCEJYyAjQjBCNEJBIiwKIlshCWMgI0IzQjZCOSIsCiJ9IQljICNBREFEQjAiLAoifCEJYyAjQUNBQkFGIiwKIjEhCWMgIzNGM0ZBNyIsCiIyIQljICM4MTdEQUEiLAoiMyEJYyAjRTFBQ0FEIiwKIjQhCWMgI0M4QUNBQyIsCiI1IQljICNGN0Y3RjciLAoiNiEJYyAjRjlGN0Y3IiwKIjchCWMgI0RBQzhDRiIsCiI4IQljICM1MDQzQUEiLAoiOSEJYyAjRTNBMEFCIiwKIjAhCWMgI0M1OTBBQSIsCiJhIQljICMyMjE5QTMiLAoiYiEJYyAjMDMwMkExIiwKImMhCWMgIzc2NjRBOCIsCiJkIQljICM4MTgzQUYiLAoiZSEJYyAjQURBRkIxIiwKImYhCWMgI0FCQUJBRiIsCiJnIQljICNBQkFCQjAiLAoiaCEJYyAjQUJBQ0IwIiwKImkhCWMgI0FDQUNBRiIsCiJqIQljICMzNjM2QTQiLAoiayEJYyAjNkU2QkE5IiwKImwhCWMgI0QwQURBRCIsCiJtIQljICNERkFDQUQiLAoibiEJYyAjQzVBQ0FCIiwKIm8hCWMgI0FEQUJBQiIsCiJwIQljICM1RjVGQTciLAoicSEJYyAjMTAxMEEyIiwKInIhCWMgIzE1MTVBOSIsCiJzIQljICNDOUM5RUMiLAoidCEJYyAjMzUzM0IxIiwKInUhCWMgIzMwMkJBNCIsCiJ2IQljICM5NzgxQUYiLAoidyEJYyAjQzFBOEJDIiwKInghCWMgI0QzQjhCNiIsCiJ5IQljICNCQkEyQUUiLAoieiEJYyAjNEM0MUE2IiwKIkEhCWMgIzNFMkNBNCIsCiJCIQljICNGMEFBQUIiLAoiQyEJYyAjQ0Y5NkFBIiwKIkQhCWMgIzU5NDJBNSIsCiJFIQljICM5NzZGQTgiLAoiRiEJYyAjQjQ4NUE5IiwKIkchCWMgIzkwNkJBNyIsCiJIIQljICM0QjM4QTQiLAoiSSEJYyAjODI4NEIwIiwKIkohCWMgI0FFQjBCMiIsCiJLIQljICNBQkFEQUQiLAoiTCEJYyAjQTJBNUIwIiwKIk0hCWMgIzlFQTFCMCIsCiJOIQljICNBRUFDQjciLAoiTyEJYyAjMUMxQ0E2IiwKIlAhCWMgIzcyNkJBOCIsCiJRIQljICNEQ0FEQUQiLAoiUiEJYyAjQzNBQkFCIiwKIlMhCWMgIzA3MDdBNCIsCiJUIQljICNCOEI4RTUiLAoiVSEJYyAjNjg2OEM4IiwKIlYhCWMgIzZBNTlBQyIsCiJXIQljICNDQ0FGQUQiLAoiWCEJYyAjQzNCMkIxIiwKIlkhCWMgI0M5QjdCOSIsCiJaIQljICNERUM1QjkiLAoiYCEJYyAjRjNENUI1IiwKIiB+CWMgI0Y1RDZCMCIsCiIufgljICM3RDY5QUIiLAoiK34JYyAjOEE2MkE3IiwKIkB+CWMgIzNGMkVBNCIsCiIjfgljICMwNzA1QTEiLAoiJH4JYyAjQUY4MEE4IiwKIiV+CWMgIzk3NzNBQSIsCiImfgljICM5Mjk2QzAiLAoiKn4JYyAjQzVDOEM2IiwKIj1+CWMgI0IyQjdCQiIsCiItfgljICNCOUJBQkUiLAoiO34JYyAjQjRCNUI4IiwKIj5+CWMgI0FEQjFCMyIsCiIsfgljICNCRUMxQzEiLAoiJ34JYyAjQzRDN0NBIiwKIil+CWMgI0JGQzhEMCIsCiIhfgljICNDQkNCQzQiLAoifn4JYyAjMUYxRkE1IiwKInt+CWMgIzc4NkNBOCIsCiJdfgljICNEOEFFQUUiLAoiXn4JYyAjREFBREFEIiwKIi9+CWMgI0MwQUJBQiIsCiIofgljICNBMkEyQUEiLAoiX34JYyAjOUU5RUFBIiwKIjp+CWMgIzk0OTRBOSIsCiI8fgljICM5OTk5Q0QiLAoiW34JYyAjRTRFNEY1IiwKIn1+CWMgI0Y1RjVGQiIsCiJ8fgljICMzMTMxQjMiLAoiMX4JYyAjRTlDQkNBIiwKIjJ+CWMgI0RGQjJCMSIsCiIzfgljICNDNkFEQUUiLAoiNH4JYyAjQzFCM0I2IiwKIjV+CWMgI0MxQjdCMiIsCiI2fgljICNEN0M5QjQiLAoiN34JYyAjRURENUI2IiwKIjh+CWMgI0YyQ0FCNCIsCiI5fgljICMyRjIzQTMiLAoiMH4JYyAjMTMwREExIiwKImF+CWMgI0JGODlBOSIsCiJifgljICM3OTU4QTYiLAoiY34JYyAjNUM0OEE2IiwKImR+CWMgIzkzOTNCOCIsCiJlfgljICNDQ0M5QjkiLAoiZn4JYyAjQzNEMUNEIiwKImd+CWMgI0NEQkRCOCIsCiJofgljICNDNEM3QzYiLAoiaX4JYyAjQjNCREJFIiwKImp+CWMgI0NDQ0FDOCIsCiJrfgljICNDQUMxQzAiLAoibH4JYyAjQkZDN0NEIiwKIm1+CWMgI0M4RDNDOSIsCiJufgljICMxRjIwQTciLAoib34JYyAjODI2REE5IiwKInB+CWMgI0Q3QURBRSIsCiJxfgljICNCRUFCQUIiLAoicn4JYyAjQjVCNUU0IiwKInN+CWMgI0EzQTNERCIsCiJ0fgljICM5Njk2RDkiLAoidX4JYyAjRjZFQUVBIiwKInZ+CWMgI0U0QkZCRSIsCiJ3fgljICNEOEFGQUYiLAoieH4JYyAjQzNBREIyIiwKInl+CWMgI0I1QURCMiIsCiJ6fgljICNCM0FFQUQiLAoiQX4JYyAjQzNCN0FFIiwKIkJ+CWMgI0UyQkVCMyIsCiJDfgljICM3QzVDQTYiLAoiRH4JYyAjRDM5NkFBIiwKIkV+CWMgIzc5NTdBNiIsCiJGfgljICNEOTlFQUIiLAoiR34JYyAjQjk4RUFBIiwKIkh+CWMgIzhGOTNCQSIsCiJJfgljICNEMEQ1QzYiLAoiSn4JYyAjQjZDM0M1IiwKIkt+CWMgI0JGQzJDNiIsCiJMfgljICNDOEMzQzIiLAoiTX4JYyAjQjlCOUJEIiwKIk5+CWMgI0FFQjFCNiIsCiJPfgljICNCQUJBQkUiLAoiUH4JYyAjQjlCNkJBIiwKIlF+CWMgI0JCQkNCRiIsCiJSfgljICNDOUM2QkQiLAoiU34JYyAjMjAxRkE2IiwKIlR+CWMgIzg4NkRBQSIsCiJVfgljICNEQ0FCQUQiLAoiVn4JYyAjN0I3QkNGIiwKIld+CWMgI0FBNzlBOCIsCiJYfgljICM1NzNGQTUiLAoiWX4JYyAjMDUwNEExIiwKIlp+CWMgIzk2OTJCQiIsCiJgfgljICNDQkNDQzAiLAoiIHsJYyAjQzVDQUM5IiwKIi57CWMgI0IzQjdCQSIsCiIrewljICNBQkFDQUYiLAoiQHsJYyAjQUJBQkIxIiwKIiN7CWMgI0FGQjBCNyIsCiIkewljICNDOEM0QzUiLAoiJXsJYyAjQ0VCQkI2IiwKIiZ7CWMgIzIyMURBNCIsCiIqewljICM4QTZCQTkiLAoiPXsJYyAjREFBQ0FEIiwKIi17CWMgI0QzQUVBRiIsCiI7ewljICNBQUFBQjMiLAoiPnsJYyAjMEMwQ0E2IiwKIix7CWMgIzBFMEVBNiIsCiInewljICMxMzExQTUiLAoiKXsJYyAjMTIwRUEyIiwKIiF7CWMgIzFGMTZBMiIsCiJ+ewljICMzQTI5QTMiLAoie3sJYyAjRDc5OUFBIiwKIl17CWMgIzUxM0FBNCIsCiJeewljICNEOUE0QUMiLAoiL3sJYyAjOUU4N0IwIiwKIih7CWMgI0Q4QzJCQiIsCiJfewljICNEMUM4QzAiLAoiOnsJYyAjQkVCQUI5IiwKIjx7CWMgI0FGQUJBRSIsCiJbewljICNBRkFCQUMiLAoifXsJYyAjQUVBQkFDIiwKInx7CWMgI0IwQUJBQyIsCiIxewljICNCNEFDQUUiLAoiMnsJYyAjQzFCNEI2IiwKIjN7CWMgI0QwQjNCMiIsCiI0ewljICMyMTFCQTMiLAoiNXsJYyAjODk2QkE4IiwKIjZ7CWMgI0Q4QUNBRCIsCiI3ewljICNEMUFFQUYiLAoiOHsJYyAjQjdBQkFCIiwKIjl7CWMgI0NBQ0FENCIsCiIwewljICNCNkI2RTQiLAoiYXsJYyAjRkNGOEY4IiwKImJ7CWMgI0VFQ0RDRCIsCiJjewljICNGMUFEQUYiLAoiZHsJYyAjRjBBREFFIiwKImV7CWMgIzY5NENBNiIsCiJmewljICNEMzlCQUEiLAoiZ3sJYyAjQUY4M0E5IiwKImh7CWMgI0E0N0ZBOSIsCiJpewljICNEM0FFQjAiLAoiansJYyAjQzFBREFFIiwKImt7CWMgI0MxQURBRCIsCiJsewljICNDM0FEQUQiLAoibXsJYyAjQ0FBQ0FEIiwKIm57CWMgI0QxQURBRSIsCiJvewljICNEN0FDQUMiLAoicHsJYyAjMjIxQkEzIiwKInF7CWMgIzg4NkJBNyIsCiJyewljICNENkFDQUMiLAoic3sJYyAjQ0NBRUIwIiwKInR7CWMgI0I0QUJBQiIsCiJ1ewljICNDN0M1QzYiLAoidnsJYyAjRUJFQkY4IiwKInd7CWMgI0Q1RDVGMCIsCiJ4ewljICNGQUYzRjMiLAoieXsJYyAjRUVDRENFIiwKInp7CWMgI0VEQjBCMiIsCiJBewljICNGNEFDQUIiLAoiQnsJYyAjRjJBQ0FDIiwKIkN7CWMgI0E1NzdBOCIsCiJEewljICM2QTRFQTYiLAoiRXsJYyAjM0UyRkE0IiwKIkZ7CWMgI0E1N0VBOSIsCiJHewljICNEN0FEQUQiLAoiSHsJYyAjRDhBQ0FDIiwKIkl7CWMgIzg4NkJBOCIsCiJKewljICNEN0FCQUMiLAoiS3sJYyAjQzhBREFFIiwKIkx7CWMgI0MwQzFDMCIsCiJNewljICM5QjlCREIiLAoiTnsJYyAjMzczN0I1IiwKIk97CWMgI0VGRDZENiIsCiJQewljICNFRUI0QjMiLAoiUXsJYyAjQzU4Q0FBIiwKIlJ7CWMgIzZGNTBBNiIsCiJTewljICM1QzQyQTUiLAoiVHsJYyAjRTVBNEFCIiwKIlV7CWMgI0U5QThBQiIsCiJWewljICMxRTE1QTIiLAoiV3sJYyAjOTE2QkE4IiwKIlh7CWMgI0UxQThBQiIsCiJZewljICM2ODRFQTUiLAoiWnsJYyAjODY2QkE4IiwKImB7CWMgI0MwQUNBQyIsCiIgXQljICMwRjBGQTciLAoiLl0JYyAjNTg1OEMxIiwKIitdCWMgIzc3NzdDRCIsCiJAXQljICM1MDUwQkUiLAoiI10JYyAjREJEQkYyIiwKIiRdCWMgI0NGQ0ZFRCIsCiIlXQljICNBNUE1REUiLAoiJl0JYyAjRDZENkYwIiwKIipdCWMgI0FFQURFMSIsCiI9XQljICM3RDcxQzAiLAoiLV0JYyAjMzAyNUE1IiwKIjtdCWMgI0QxOTZBQSIsCiI+XQljICNBQTdDQTgiLAoiLF0JYyAjMzgyQUE0IiwKIiddCWMgIzc2NTdBNyIsCiIpXQljICM4RjZCQTgiLAoiIV0JYyAjMjYxQ0EyIiwKIn5dCWMgI0RGQUFBQiIsCiJ7XQljICNEN0E3QUIiLAoiXV0JYyAjQ0RBMUFBIiwKIl5dCWMgI0NDQTFBQSIsCiIvXQljICNDQkExQUEiLAoiKF0JYyAjMUYxOUEzIiwKIl9dCWMgIzdENjVBNyIsCiI6XQljICNDOEExQUEiLAoiPF0JYyAjQzNBMkFDIiwKIltdCWMgI0I4QUFBQyIsCiJ9XQljICM2QzZDQzkiLAoifF0JYyAjQzJDMkU5IiwKIjFdCWMgIzQ0MzJBNSIsCiIyXQljICNCRDg3QTkiLAoiM10JYyAjRURBQkFEIiwKIjRdCWMgIzk5NzBBOCIsCiI1XQljICMxRTE3QTMiLAoiNl0JYyAjQTA3OUE4IiwKIjddCWMgI0M4OUJBQSIsCiI4XQljICM5MjczQTgiLAoiOV0JYyAjMEUwREExIiwKIjBdCWMgI0E5QTlCMSIsCiJhXQljICNEMEQwRDEiLAoiYl0JYyAjREZERkY0IiwKImNdCWMgIzE3MTdBOSIsCiJkXQljICNGM0YzRkIiLAoiZV0JYyAjODg4OEQzIiwKImZdCWMgIzNCM0JCNyIsCiJnXQljICMxMzEzQTgiLAoiaF0JYyAjNzg3OENEIiwKImldCWMgI0RFRDhFQSIsCiJqXQljICNFRUQyRDEiLAoia10JYyAjRTlBRkFGIiwKImxdCWMgI0NFOThBQSIsCiJtXQljICM2NTRBQTUiLAoibl0JYyAjMzMyNkE0IiwKIm9dCWMgIzdGNUZBNyIsCiJwXQljICNEQ0E1QUIiLAoicV0JYyAjRDVBMUFBIiwKInJdCWMgIzBCMDhBMSIsCiJzXQljICMxODEyQTIiLAoidF0JYyAjQ0I5REFBIiwKInVdCWMgIzk5NzhBOCIsCiJ2XQljICMxNzE2QTMiLAoid10JYyAjRDhEOERGIiwKInhdCWMgIzYzNjNDNiIsCiJ5XQljICM5MDkwRDYiLAoiel0JYyAjRUVFRUY5IiwKIkFdCWMgIzA2MDZDRiIsCiJCXQljICNGOEVFRUUiLAoiQ10JYyAjRUNDOEM4IiwKIkRdCWMgI0U3QjBCMiIsCiJFXQljICNEQUFFQUQiLAoiRl0JYyAjRENBRkFEIiwKIkddCWMgI0REQUZBRCIsCiJIXQljICNFMEFFQUMiLAoiSV0JYyAjRDRBQ0FCIiwKIkpdCWMgI0NGQUNBRSIsCiJLXQljICNDQ0IyQjMiLAoiTF0JYyAjREVEMUQxIiwKIk1dCWMgI0Y4RjdGNyIsCiJOXQljICMzMzMzRDciLAoiT10JYyAjRThFOEZCIiwKIlBdCWMgI0ZCRjdGOCIsCiJRXQljICNFRkQ4RDkiLAoiUl0JYyAjREFCN0I4IiwKIlNdCWMgI0JGQURBQyIsCiJUXQljICNCRUFDQUMiLAoiVV0JYyAjQzJBREFEIiwKIlZdCWMgI0M1QURBRCIsCiJXXQljICNENkFFQUYiLAoiWF0JYyAjRERBRUFGIiwKIlldCWMgI0UwQUNBRSIsCiJaXQljICNEMUFDQUIiLAoiYF0JYyAjQ0VBQ0FFIiwKIiBeCWMgI0QyQkJCRCIsCiIuXgljICNFN0RFREUiLAoiK14JYyAjRkNGQkZCIiwKIkBeCWMgIzFFMUVENCIsCiIjXgljICNCOEI4RjEiLAoiJF4JYyAjRkFGOUY5IiwKIiVeCWMgI0Q3RDBDRiIsCiImXgljICNDMkJGQkUiLAoiKl4JYyAjQjZCOEI4IiwKIj1eCWMgI0I5QjhCQiIsCiItXgljICNCQkI0QjciLAoiO14JYyAjQjhCMkIyIiwKIj5eCWMgI0I0QjFBRiIsCiIsXgljICNCNEIwQUYiLAoiJ14JYyAjQkJBRUFFIiwKIileCWMgI0NEQjFCMSIsCiIhXgljICNENkFGQjAiLAoifl4JYyAjRENBQ0FEIiwKInteCWMgI0QxQURBRCIsCiJdXgljICNDRUFFQUUiLAoiXl4JYyAjRDJCOUI5IiwKIi9eCWMgI0UxRDRENCIsCiIoXgljICNGQkZBRkEiLAoiX14JYyAjRjZGNkZFIiwKIjpeCWMgIzA4MDhDRiIsCiI8XgljICM2NzY3RTIiLAoiW14JYyAjREZFMkUyIiwKIn1eCWMgI0M5Q0JDRSIsCiJ8XgljICNCOUJBQkQiLAoiMV4JYyAjQjlCREJFIiwKIjJeCWMgI0I0QjhCOCIsCiIzXgljICNCOEI5QjkiLAoiNF4JYyAjQkNCQkJDIiwKIjVeCWMgI0JFQjhCQSIsCiI2XgljICNCQ0FGQjIiLAoiN14JYyAjQkJBQkFDIiwKIjheCWMgI0NGQURBRCIsCiI5XgljICNDQ0IxQjEiLAoiMF4JYyAjRDZDM0M0IiwKImFeCWMgI0YxRUJFQiIsCiJiXgljICNBNEE0RUUiLAoiY14JYyAjMjAyMEQ1IiwKImReCWMgI0ZFRkVGRiIsCiJlXgljICNGMEYyRjIiLAoiZl4JYyAjRDREOEQ5IiwKImdeCWMgI0NBQ0NDRSIsCiJoXgljICNDQUM5Q0IiLAoiaV4JYyAjQzNDN0M4IiwKImpeCWMgI0JGQzNDNiIsCiJrXgljICNCQUJEQzAiLAoibF4JYyAjQjRCM0I1IiwKIm1eCWMgI0I5QUVBRiIsCiJuXgljICNDM0FDQUMiLAoib14JYyAjRDFBRUFFIiwKInBeCWMgI0Q1QUVBRSIsCiJxXgljICNENUFGQjAiLAoicl4JYyAjRDZCNUI1IiwKInNeCWMgI0Q5QkNCQyIsCiJ0XgljICNEQkMxQzEiLAoidV4JYyAjREZDOUM5IiwKInZeCWMgI0U5RENEQyIsCiJ3XgljICNGNEVERUQiLAoieF4JYyAjRkJGOUY5IiwKInleCWMgIzUwNTBERSIsCiJ6XgljICM5MjkyRUEiLAoiQV4JYyAjRUZGMEYwIiwKIkJeCWMgI0U3RThFOCIsCiJDXgljICNFMUUwRTEiLAoiRF4JYyAjREVEQ0REIiwKIkVeCWMgI0UxREJEQyIsCiJGXgljICNFMUQ5RDkiLAoiR14JYyAjRTNENkQ3IiwKIkheCWMgI0UzRDNENCIsCiJJXgljICNFN0Q3RDciLAoiSl4JYyAjRThEQ0RDIiwKIkteCWMgI0VERTVFNiIsCiJMXgljICNGNUYwRjAiLAoiTV4JYyAjRDFEMUY3IiwKIk5eCWMgIzA1MDVDRiIsCiJPXgljICMxOTE5RDMiLAoiUF4JYyAjRUZFRkZDIiwKIlFeCWMgIzNFM0VEQSIsCiJSXgljICM3MzczRTQiLAoiU14JYyAjNTU1NURFIiwKIlReCWMgIzhBOEFFOCIsCiJVXgljICMwMjAyQ0YiLAoiVl4JYyAjMzYzNkQ4IiwKIldeCWMgI0U1RTVGQiIsCiJYXgljICNGOEY4RkUiLAoiWV4JYyAjRjBGMEZDIiwKIlpeCWMgIzNDM0NEOSIsCiJgXgljICM3QTdBRTUiLAoiIC8JYyAjOUI5QkVCIiwKIi4vCWMgI0I2QjZGMiIsCiIrLwljICNCQUJBRjEiLAoiQC8JYyAjQTFBMUVEIiwKIiMvCWMgIzg1ODVFNyIsCiIkLwljICM0QjRCREQiLAoiJS8JYyAjMEUwRUQxIiwKIiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAiLAoiICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICIsCiIgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAuICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIiwKIiAgICAgICAgICAgICAgICAgICAgICAgIC4gKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgICAgICAgICAgICAgICAgICAgICAgICAgICAiLAoiICAgICAgICAgICAgICAgICAgICAuICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgICAgICAgICAgICAgICAgICAgICAgICIsCiIgICAgICAgICAgICAgICAgICArICsgKyArIEAgIyAkICUgJiAqICogKiAqICogKiAqICogKiAqICogKiAqICogKiAqICogKiAqICogKiAqICogKiAqICogKiAqICogKiAqICogKiAqICogKiAqICogKiAqICogKiAqICogKiAqICogKiAqICogKiAqICogKiAqICogKiAqICogPSAtIDsgPiAsICsgKyArICsgICAgICAgICAgICAgICAgICAgIiwKIiAgICAgICAgICAgICAgJyArICsgKyApICEgfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IHsgXSArICsgKyAnICAgICAgICAgICAgICAgICAiLAoiICAgICAgICAgICAgICArICsgKyBeIC8gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+ICggXyArICsgKyAgICAgICAgICAgICAgICIsCiIgICAgICAgICAgICArICsgKyA6IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiA8IFsgKyArICAgICAgICAgICAgICAgIiwKIiAgICAgICAgICArICsgKyB9IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfCAsICsgKyAgICAgICAgICAgICAiLAoiICAgICAgICAgICsgKyAxIDIgfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IF4gKyArIC4gICAgICAgICAgICIsCiIgICAgICAgICsgKyArIDMgfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiA0IDUgNiA3IDggOSAwIGEgYiBjIGQgZSBmIGcgaCBpIGogayBsIG0gNCBuIH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gbyBwICsgKyAgICAgICAgICAgIiwKIiAgICAgICAgKyArIHEgfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gciBzIHQgdSB2IHcgeCB5IHogQSBCIEMgRCBFIEYgRyBIIEkgSiBLIEwgTSBOIE8gUCByIH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IFEgKyArICAgICAgICAgICAiLAoiICAgICAgICArICsgUiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gbiBTIFQgVSBWIFcgWCBZIFogYCAgLi4uKy5ALiMuJC4lLiYuKi5IID0uLS47Lj4uLC4nLikuIS5+LnsuXS5yIH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gXi4rICsgLiAgICAgICAgICIsCiIgICAgICArICsgKyAvLn4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gNCAoLl8uOi48LlsufS58LjEuMi4zLjQuNS42LjcuOC45LjAuYS5iLmMuZC5lLmYuZy5oLmkuai5rLmwubS5uLm8ucC5+IH4gfiB+IH4gcS5yLnMudC51LnQucy52LncueC55LnouQS5+IH4gfiB+IH4gfiB+IH4gfiBCLkMuKyArICAgICAgICAgIiwKIiAgICAgICsgKyArIEQufiB+IH4gRS5GLnUuRy5HLkcuRy5HLkcuRy5HLkcuRy5HLkcuRy5HLkcuSC5JLjEuMS4xLkouSy5MLksuTS5OLk8uUC5HLkcuRy5HLkcuRy5HLkcuUS5ILkcuUi5TLlQuVC5VLlYuVy5YLlkuNCB+IFouYC5HLkcuRy5HLkcuRy5HLkcuRy4gKy4rRy4rK34gfiB+IH4gfiB+IH4gfiB+IH4gQCsrICsgICAgICAgICAiLAoiICAgICAgKyArICMrfiB+IH4gfiAkK0cuRy5HLkcuRy5HLkcuRy5HLkcuRy5HLkcuRy5HLkcuRy5HLkcuJSsmKyYrKitMLj0rPSstKzsrRy5HLkcuPissKycrKSssK0cuRy5HLkcuRy4hK34reytdK3srXisvKygrXys6KzwrWytHLkcuRy59K3wrMStzLjIrRy5HLkcuRy5HLjMrfiB+IH4gfiB+IH4gfiB+IH4gfiA0KysgKyAgICAgICAgICIsCiIgICAgICArICsgNSt+IH4gfiB+IDYrNys4KzkrRy5HLkcuMCswKzArMCswK2ErYitjK2QrZStmK0cuRy5nK2graStqK2srbCtsK20rRy5HLkcubitvK3ArcCtxK3Ircyt0K0cuRy5HLkgudSt2K3Yrdit3K3greSt6K0ErQitHLkcuQytEKy8gfiB+IH4gRStGK1srRy5HLkcuMit+IH4gfiB+IH4gfiB+IH4gfiB+IEcrKyArICAgICAgICAgIiwKIiAgICAgICsgKyA1K34gfiB+IH4gfiB+IH4gSCtHLkcuWyt+IH4gfiA0IEkrSitLK0wrTStOK08rRy5HLlArUStMLlIrbCtTK1QrVStHLkcuVitXK1grWStxK3IrWitTLmArIEBHLkcuRy4uQCtAQEAjQCRAJUAmQCpAPUBHLkcuLUA7QH4gfiB+IH4gfiB+IH4gPkBbK0cuRy4sQH4gfiB+IH4gfiB+IH4gfiB+IH4gRysrICsgICAgICAgICAiLAoiICAgICAgKyArIDUrfiB+IH4gfiB+IH4gfiBIK0cuRy5bK34gfiAnQClAIUB+QHtATC5dQGsrXkBHLkcuL0AoQFIrUytfQDpAPEBHLkcuW0B9QHxAMUAyQDNAM0A0QFQuNUA2QEIrRy5HLjdAOEA5QDBAYUBiQEBAY0BHLkcuZEBlQGZAfiB+IH4gfiB+IH4gfiB+IGdARy5HLixAfiB+IH4gfiB+IH4gfiB+IH4gfiBHKysgKyAgICAgICAgICIsCiIgICAgICArICsgNSt+IH4gfiB+IH4gfiB+IEgrRy5HLlsrfiBuIGhAaUBqQGtAbEBTK2wrUiteQEcuRy4vQFIrbUBfQDpAOkBuQEcuRy5vQFgrWCtxK3BAcitxQDRAckBzQHRAdUBHLnZAd0B4QHlAekBBQEFAQkBHLkcuRy5DQERARUBGQH4gfiB+IH4gfiB+IH4gR0AsQEcuSEB+IH4gfiB+IH4gfiB+IH4gfiB+IEcrKyArICAgICAgICAgIiwKIiAgICAgICsgKyA1K34gfiB+IH4gfiB+IH4gSCtHLkcuWysnQElASkBLQExATUBOQE9AUytsK1BARy5HLlFAUkBTQFRAX0BfQFVARy5HLlZAWStZKzNAcUAzQFdAWEA1QFlAWkBgQCAjLiMrI0AjIyMkIyUjJiMqI0cuRy49Iy0jOyM+IywjJyN+IH4gfiB+IH4gfiB+ICkjISN+I34gfiB+IH4gfiB+IH4gfiB+IH4gRysrICsgICAgICAgICAiLAoiICAgICAgKyArIDUrfiB+IH4gfiB+IH4gfiBIK0cuRy5bK3sjXSNeIy8jKCNfI0cuOiNTQFNAPEBHLjwjUkBbI30jfCMxIzpAMiNHLkcuZiszI3IrcUBxQHFANCNYQDVAWUA1IzYjNyMrIzgjOSMwI2EjYiNjI0IrRy5HLmQjZSNmI2cjaCNpI2ojfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiBHKysgKyAgICAgICAgICIsCiIgICAgICArICsgNSt+IH4gfiB+IH4gfiB+IEgrRy5HLlsrayNsI20jbiNvI3AjRy5xI1RAU0BsK3IjcyN0I1RAdSN1I30jdiN3I0cuRy5HLngjeSNaK1orUy40I1hAeiNaQEEjQiNDI0QjOSMwI2EjRSNGI0cjRy5HLkcuSCNJI0ojSyNMI00jTiNPI34gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IEcrKyArICAgICAgICAgIiwKIiAgICAgICsgKyA1K34gfiB+IH4gfiB+IH4gSCtHLkcuLitQI1EjUiNTI1QjdkBHLlUjViM6QFMrdCNTK19AUytXI24jWCNZI1ojYCNHLkcuRy5CKyAkLiQrJEAkIyRZQCQkJSQ3IyYkKiQwI2IjPSQtJDskPiRHLkcuRy4sJCckKSQhJH4keyRdJF4kfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gRysrICsgICAgICAgICAiLAoiICAgICAgKyArIDUrfiB+IH4gfiB+IH4gfiBIK0cuRy5HLkcuRy5HLkcuRy5HLkcuVSNXK1YjX0BTQFNALyQoJF8kOiQ8JFskfSR8JDEkRy5HLkcuRy5HLkcuRy5HLjIkMyQ0JDUkNiQ5IzckOCRGIzkkMCRhJEcuRy5iJGMkZCRlJGYkZyRoJF0kaCRpJH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiBHKysgKyAgICAgICAgICIsCiIgICAgICArICsgNSt+IH4gfiB+IH4gfiB+IEgrRy5HLkcuRy5HLkcuRy5HLkcuRy5VI1crViNqJHwjayRsJG0kbiRvJHAkcSRyJHMkdCR1JHYkRy5HLkcuRy5HLkcuRy5HLkcudyR4JHkkeiQ9JEEkQiRDJEQkRy5HLkUkRiRHJEgkSSRKJEskcyRMJE0kTiR+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IEcrKyArICAgICAgICAgIiwKIiAgICAgICsgKyA1K34gfiB+IH4gfiB+IH4gTyRHLkcuUCRRJFIkUyRUJFUkPCNHLlYkVytWI2okVyRYJFkkWiRgJCAlLiUrJUAlSyQjJSQlJSUmJSolPSUtJUIrRy5HLkcuRy5HLkcuOyU+JSwlJyUpJSElfiVHLkcueyVdJV4lLyUoJV8lOiU8JVslfSV8JTElMSUxJTElMSUxJTElMiV+IH4gfiB+IH4gfiB+IH4gRysrICsgICAgICAgICAiLAoiICAgICAgKyArIDUrfiB+IH4gfiB+IH4gciAzJUcuRy5QJDQlNSU2JTclOCU5JUcubkB8QFcrdiMwJWElYiVjJWQlZSVmJWgkZyVoJWklcyRjJTYlaiVrJWwlbSVuJW8lcCVCK0cuRy5HLnElciU7JHMldCV1JUcuRy5HLnYldyU5IHgleSV6JU0kQSVHLkcuRy5HLkcuRy5HLkcuRy5CJX4gfiB+IH4gfiB+IH4gfiBHKysgKyAgICAgICAgICIsCiIgICAgICArICsgNSt+IH4gfiB+IH4gfiBDJUQlRy5HLlAkRSVJIEYlRyVIJUklSiVLJVgrVytMJU0lTiU9Lk8lUCVRJVAlUCVSJWgjaCNTJWclVCVVJVYlVyVYJVklWiVgJSAmLiZHLkcuKyZyJTskQCYjJiQmRy5HLkcuJSYmJiomPSYtJjsmPiYsJlAkRy5HLkcuRy5HLkcuRy4nJikmfiB+IH4gfiB+IH4gfiB+IEcrKyArICAgICAgICAgIiwKIiAgICAgICsgKyA1K34gfiB+IH4gfiB+ICEmfiZHLkcuUCR7Jl0mXiYvJigmXyZxKzJAMUA6JjwmWyZ9JnwmMSYyJmMlUCVoJTMmNCY1JjYmaCU3JjgmOSZIJTAmLytZQGEmYiZjJkcuRy5kJnIlOSRAJiMmZSZHLkcuRy5mJmcmaCZpJmomayZsJmYkKyVdJG0mbiZ+IG8mRy5HLnAmfiB+IH4gfiB+IH4gfiB+IH4gRysrICsgICAgICAgICAiLAoiICAgICAgKyArIDUrfiB+IH4gfiB+IG4gcSZyJkcuRy4+K3MmNy50JnUmdiZ2JjNAM0AxQG8jdyZHLkcueCZ5JnomQSZSJUImUCRHLkMmRCZoI0UmUSRGJmxARyZYJVlAJCRiJkgmRy5HLkIrSSYnJUomSyZMJk0mRy5HLk4mTyZQJlEmUiZTJlQmVSZnJHMkViZXJn4gbyZHLkcucCZ+IH4gfiB+IH4gfiB+IH4gfiBHKysgKyAgICAgICAgICIsCiIgICAgICArICsgNSt+IH4gfiB+IH4gWCZZJlomRy5HLj4rYCY7LiAqLipULnFAM0AzQCsqQCojKkcuRy4kKiUqJioqKj0qLSpHLkcuRy47Kj4qLConKikqISp+KjVAeyokJF0qXipHLkcuLyooKl8qMCQ6KjwqWypHLkcuRy59KnwqUCYxKjIqMyo0KjUqNipdJDUqJ0BvJkcuRy5wJn4gfiB+IH4gfiBuIDcqfiB+IEcrKyArICAgICAgICAgIiwKIiAgICAgICsgKyA1K34gfiB+IH4gciA4KjkqMCpHLkcuUCRhKmIqYyo0I2QqZCpxQDNAai5lKmYqRy5HLmcqfCZoKmkqaiprKkcuRy5HLmwqbSosKm4qbypwKlUuNUB7KiQkYCVxKkcuRy5yKnMqdCp1KnYqdyp4KnkqRy5HLlsreipBKkIqQypEKkUqZiRGKl0kYyU3Km8mRy5HLnAmfiB+IH4gfiBHKkgqSSpOJH4gRysrICsgICAgICAgICAiLAoiICAgICAgKyArIDUrfiB+IH4gfiBKKksqTCpNKkcuRy5QJE4qTypYJVhAUCpkKldAM0AzQFEqUipHLkcuUypUKkYqVSpWKlcqRy5HLkcuQitYKj0uWSpaKmAqXyYvK3ojWUAgPS4mRy5HLi49Kz1yJTkkSiZAPSM9JD0uK0cuRy4lPSY9Kj09PS09Oz0+PUskYyVMJCw9byZHLkcucCZ+IH4gfiAnPUwkZyQpPSE9fiBHKysgKyAgICAgICAgICIsCiIgICAgICArICsgNSt+IH4gfiB+IH49ez1dPV49Ry5HLkIrLz0oPV89Xz1fPTo9Oj06PTo9PD1bPUcuRy59PXw9MT0yPVQqMz1HLkcuRy5HLjQ9NT02PTc9OD12JlUuOT0wPWE9Ry5HLmI9Yz1kPWU9QSRmPWc9aD1pPWo9Ry5HLkcuRy5rPWw9bT07PVUmbSZuPW89cD1CK0cuRy5wJn4gJ0BxPTYqTCRWJmckTSR+IEcrKyArICAgICAgICAgIiwKIiAgICAgICsgKyA1K34gfiB+IH4gcj1HLkcuRy5HLkcuRy5HLkcuRy5HLkcuRy5HLkcuRy5HLkcuRy5HLnM9dD1nJUwkaCMsKnU9Ry5HLkcuRy5HLkgudj13PXg9eT16PUcuRy5HLkE9Qj1DPUQ9RT1GIzskcyVLJkY9Rz1IPXslRy5HLkcuRy5JPUo9Sz1MPUIrRy5HLkcuRy5HLk09KT1nJEwkcyRzJF0kcyROPX4gRysrICsgICAgICAgICAiLAoiICAgICAgKyArIDUrfiB+IH4gfiBPPVA9Ry5HLkcuRy5HLkcuRy5HLkcuRy5HLkcuRy5HLkcuRy5ILlE9TCRSPVM9aCNoI1Q9VT1HLkcuVj1QJEcuRy5HLkcuRy5HLkcuRy5HLlc9WD1ZPVo9YD0gLS4tdCorLUAtIy0kLSUtJi0qLUcuRy5HLkcuRy5HLkcuRy5HLkcuRy49LS0tOy1WJl0kPj0tPS09cT0tPTcqfiBHKysgKyAgICAgICAgICIsCiIgICAgICArICsgNSt+IH4gfiB+ID4tLC0nLSktIS0hLX4tey1dK10teyteLX4rckBULlMuM0AvLVckKC1gIF8tOi1pJUwkPC1bLUcufS18LTEtMi0zLTQtNS02LTctOSU4LTktMC1hLWItYy1kLWUtZi1nLWgtaS1qLWstbC1tLW4tby1wLXEtci1zLXQtdS12LXcteC15LT49SyQrJT4mei0nQH4gfiB+IH4gfiB+IEcrKyArICAgICAgICAgIiwKIiAgICAgICsgKyA1K34gfiB+IH4gQS1CLUMtRC1FLUYtRy1ILUgtSS1dLX4rSi1QKmQqNEA0QHFAWStLLUwtTS1OLS0uTy1QLVEtUi1TLVQtVS1tQHQjUytWLWokaiR1I1ctVy1YLVktWi1gLSA7LjsrO0A7IzskOyU7JjsqOz07Tj0tOy07XiQ7OyE9PjssO2gkXSRSJTUqNSpJKi09biB+IH4gfiB+IH4gfiB+IH4gRysrICsgICAgICAgICAiLAoiICAgICAgKyArIDUrfiB+IH4gfiAnOyk7ITt+O3s7XTteO0gtSC1JLS87KDtKLVAqNCNaK3FAV0AvLV87dSY6Ozw7Wzt9O307fDsxOzI7MztNK00rNDs1OzY7ais3O2srLStSQFUkODs5OzA7YTtiO2M7ZDtlO2Y7ZztoO2k7ajtIKms7aztIKlUmPiY1KkskaCRdJFIlbDttOydAfiB+IH4gfiB+IH4gfiB+IH4gfiBHKysgKyAgICAgICAgICIsCiIgICAgICArICsgNSt+IH4gfiB+IG47KTtvO3A7Jy1dO3E7SC1yO10rXi01QFhANCNXQHM7dDs4PXU7dTs6Jm8jdjt3O3g7S0BeI20jeTt6O0E7QjtDO0M7QztEO0U7RjtHO0g7STtKOzo7SztMO007TjtPO1A7UTtSO1M7VDtVO1Y7Vzs7LWkjOy1YO20mbDtZO10kTCReJG4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IEcrKyArICAgICAgICAgIiwKIiAgICAgICsgKyA1K34gfiB+IH4gWjtgOyg7ckA0QGQqUCpkKjJAMUAxQFcrViM6QDpAX0BTK2wrbCtSK2srNzs3O2orais2OzY7ID4uPis+QD4jPiQ+JT4kPiY+Iz4jPiM+Kj49Pi0+Oz4+Piw+Jz4pPiE+fj57Pl0+Xj4vPig+Xz46Piw7PD5tJiw7ZyQrJWMlTCRMJCc9biB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gRysrICsgICAgICAgICAiLAoiICAgICAgKyArIDUrfiB+IH4gfiBbPn0+fD58PjE+Mj4zPjQ+NT42Pjc+OCU4Pjk+MD5hPm4jYj5iPmM+Yz5kPmQ+MzszO2U+Zj5nPkY7aD5pPmo+az5sPmo+bT5uPm8+cD5xPj0+Oz5yPnM+dD51PnY+bCQpPiE+bCN3Png+eT56PkE+bDtCPmckbDtAJWgkViZMJFIlKT1uIH4gfiB+IH4gfiB+IH4gfiB+IH4gfiBHKysgKyAgICAgICAgICIsCiIgICAgICArICsgNSt+IH4gfiB+IEM+RD5FPkY+ITtHPkg+ST5KPncrSz5MPk0+eyo1QC4qdiZzOy8tOD11O1ctVy1OPk8+Tz5MQFJATSsxLlA+Oz5RPlI+PT47Pi0+Uz5UPlQ+VT5WPlc+WSNEO1g+WT4yO1o+YD4gLC4sKyxALCMsJCxAJUAlQCUrJWgkcyRdJEwkTCQlLH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IEcrKyArICAgICAgICAgIiwKIiAgICAgICsgKyA1K34gfiB+IH4gJzsmLCosPSwtLC0sfjtxO0MtOyxJLUktXS0+LD4sckBULlQuUy5aK3IrcStZK3ArcCs6JjpAUytSKywsJisnLCs+Kz4nLEM7QjsuPiksKSx8LiEsfC5+LHssXSwzO14sLywoLF8sOiw8LFssfSx8LHwsMSw2KmgkWTtjJV0kTCRSJSc9fiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gRysrICsgICAgICAgICAiLAoiICAgICAgKyArIDUrfiB+IH4gfiAyLDMsNCw9LD0sNSxELTYsRy17LXstNyxdLT4sfityQFQuZCpTLnFAMkAxQFgrWCtXK1crViNfQGwrTissLEsuNTs4LDU7OSwwLFErNjs2Ozc7KEBYI2EsYixjLGQsMD5lLGYsYDtnLGgsaSxqLE8lTyVrLGgkWTtjJVYmTCRSJUwkOz1PI34gfiB+IH4gfiB+IH4gfiB+IH4gfiBHKysgKyAgICAgICAgICIsCiIgICAgICArICsgNSt+IH4gfiB+IGwsbSxuLG8scCxxLHIsSj5eO3stIS1dLV0tPiw+LFAqUCpkKjRAM0AyQDFAWCtXK1crVytXK3wjU0BtQGwrcyx0LChAaytSKy0rdSxjPnYsdyxsQExAVyM6QDomeCx5LHosWS0wLUEsQixDLFk7RCZEJlk7cyRjJVYmXSRMJFIlXSRmJG4tRCxqI34gfiB+ICcjKT1FLHwqfiB+IEcrKyArICAgICAgICAgIiwKIiAgICAgICsgKyA1K34gfiB+IH4gRixHLEgsSSxKLEssWCBMLCQkTSxJLV0rTix7K34rSi1QKmQqNEBxQDNAMUAxQFgrVytXK1crfCNPLHQjUytqJFNAdCNQLFEsUixTLFQsVSxWLFMsVyxYLG9AWSxaLGAsICcuJysnQCcjJyQncyRzJCUnYyVjJVYmXSRMJFIlTCRWJiw7Jic7LTstVSYsPV4kbDtnJDs7KCUnQH4gRysrICsgICAgICAgICAiLAoiICAgICAgKyArIDUrfiB+IH4gfiAqJz0nLSc7Jz4nLCcnJyknISd+J00sSS0vO10rfitKLVAqZCo0QHFAM0AyQDFAWCtYK1crOkAvJE8sX0A6QDpAbyNtI3snXSdeJy8nYCBgICgnXyc6JzwnRy5HLlsnfSdzO34qfCcxJzInMyc0JzUnViZWJl0kXSRMJFAlUiVoJWckbiZ+IH4gTiReJGYkRipWJlk7SCpYO2ZAfiBHKysgKyAgICAgICAgICIsCiIgICAgICArICsgNSt+IH4gfiB+IDYnNyc4JzknMCdhJ2InYydkJ2UnZidJLTcsXS0+LEotUCpkKjRAcUAzQDJAMUAxQFgrVytnJ0xAaCc2LUcuRy5RLmknaidrJ2wnbSduJ28nbSdwJ3EncidHLkcucyc4PX0nLy1XQHQndSd2J3cndyddJF0keCdMJFAlaCVSJWglVzt+IH4gfiB+IH4gQSpeJGMlUiVdJCw7eSd+IEcrKyArICAgICAgICAgIiwKIiAgICAgICsgKyA1K34gfiB+IH4geidBJ0InQydEJ0UnRidHJ0gnSSdKJ0snSS0vO3srfit+K2QqZCpxQDNAM0AxQDFAMUAxQHArOkBMJ0cuRy5HLkcuTSdOJ08nUCdRJ1ElUSVSJ1MnVCdVJ1YnVydYJ08+dTt8QGouOD1ZJ1onOyM7I0wkTCRQJWAnICkuKVIlaCVFLH4gfiB+IH4gfiB+IH4gTj1AJWckLT0nQH4gRysrICsgICAgICAgICAiLAoiICAgICAgKyArIDUrfiB+IH4gfiBuICspQCkjKSQpJSkmKSopPSktKTspPil4KzcsPiw+LEotUCpkKjRAcUAzQDJAMUBYK1grMUA6QFcjLCknKSkpRy4hKX4peyldKV4pLykoKV8pOik8KVspfSl8KTEpdSw6JjIpWCtfOzMpNCk7I0RAUSVoJTUpRy5HLjYpUiVZOzcpfiB+IH4gfiB+IH4gfiAnQDEqbiZ+IH4gfiBHKysgKyAgICAgICAgICIsCiIgICAgICArICsgNSt+IH4gfiB+IH4gOCk5KTApYSliKWMpZCllKWYpZyloKWkpNyxLPj4sSi1QKjRANEBxQDNAMkAxQFgrVys6JjpAX0BTK20jailHLmspbCltKW4pbylwKXEpcilzKXQpdSl2KXcpYCB4KXUjZydYK2okeSk0KTI9REAyPWgleilHLkcuQSlSJTQqbiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IEcrKyArICAgICAgICAgIiwKIiAgICAgICsgKyA1K34gfiB+IH4gfiB+IEIpQylEKUUpRilHKUgpSSlKKUspTCldK3srPixKLVAqNEBTLnFAM0AyQDFAWCs6JnUjaiRXI2xATSlqKUcuTilPKVApUSlSKVMpVClVKVYpVylYKVkpWilgKSAhMSNTKzpAU0AuITYlUiUyPVIlaCUrIUcuRy5BKWMlbi1+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gRysrICsgICAgICAgICAiLAoiICAgICAgKyArIDUrfiB+IH4gfiB+IH4gQCEjISQhJSEmISohPSEtITshPiEsIXojXis+LEotUCo0QFMucUAzQCchKSEhIX4hNi17IV0hXiEvIWEuRy4oIV8hOiE8IVshfSF8ITEhQitHLkcuRy50LTIhIz50I3QjUyszITQhNiVSJTI9UiVoJXopRy5HLkEpJic1IX4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiBHKysgKyAgICAgICAgICIsCiIgICAgICArICsgNSt+IH4gfiB+IH4gfiA2ITchOCFHLkcuRy5HLkcuRy5HLkcub0A5IT4sSi1QKmQqNEBxQDAhYSFHLkcuRy5HLkcuRy5iIWMhYS5HLmQhZSFmIWchaCFpIVInaiFHLkcuRy5HLkcuayFsIW1AbUBsK20hbiFvIVIlcCFxIUIrRy5HLkcuRy5HLkcuRy5HLkcuciFzIX4gfiB+IH4gfiB+IH4gfiB+IEcrKyArICAgICAgICAgIiwKIiAgICAgICsgKyA1K34gfiB+IH4gfiB+IDtAdCFHLkcudSF2IXcheCF5IXohRy5HLkEhQiFKLVAqZCo0QEMhfiFHLkIrRCFFIUYhRyFIIUcuRy5CK0cuSSFKIWYhZiFmIWYhUSdLIUwhTSFOIU8hRy5QIU4uUitSQD0rUSFSIW8hUiVbLUIrRy5HLkcuRy5HLkcuRy5HLkcuRy5TIVQhfiB+IH4gfiB+IH4gfiB+IH4gRysrICsgICAgICAgICAiLAoiICAgICAgKyArIDUrfiB+IH4gfiB+IH4gVSFHLkcuViFXIVghWSFaIWAhIH4ufkcuRy4rfj4sUCpkKjRAQH5HLiN+JH4xQFcrViNqJGokJX5CK0cuRy4mfip+PX4tfjolO34+fix+J34pfiF+fn5HLnt+XX5rK2srTC5efi9+byFSJTI9KH5ffjp+Ry5HLjx+W35bflt+W35bfngufiB+IH4gfiB+IH4gfiB+IH4gfiBHKysgKyAgICAgICAgICIsCiIgICAgICArICsgNSt+IH4gfiB+IH4gfX4lPUcufH4xfjJ+M340fjV+Nn43fjh+OX5HLjB+PixQKmQqYX5HLkcuYn4yQDFAcCs6JjomVi19I2N+Ry5HLmR+ZX5mflYpZ35ofml+an5rfmx+bX5ufkcub34zOzc7NzswLHB+cX5vIVIlUiVSJVIleilHLkcucn5+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IEcrKyArICAgICAgICAgIiwKIiAgICAgICsgKyA1K34gfiB+IH4gfiBzfkcuRy50fnV+dn53fnh+eX56fkF+Qn5DfkcuRy5EflAqZCpFfkcuRy5GfjJAMUBYK1grVytWIzpAR35HLkcuSH5Jfkp+S35Mfk1+Tn5PflB+UX5SflN+Ry5UfksuaitRK1V+WiNqJW8hUiVSJVIlUiV6KUcuRy5yfn4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gRysrICsgICAgICAgICAiLAoiICAgICAgKyArIDUrfiB+IH4gfiB+IFZ+Ry5HLkcuRy5HLkcuRy5HLkcuRy5HLkcuRy5HLld+UCpkKlh+Ry5bQDNAMkAxQFgrViNnJzpAdSNPQFl+Ry5afmB+IHsueyt7ZiFmIUB7I3skeyV7JntHLip7Sy5KLjU7PXste00taCVSJVIlUiVSJTt7Ry5HLnJ+fiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiBHKysgKyAgICAgICAgICIsCiIgICAgICArICsgNSt+IH4gfiB+IH4gRitHLkcuPnt9Kyx7J3speyl7KXspe0E9MH4he357e3tQKmQqXXtHLltAMkAyQDFAMik6QC8jVEB9I157SC5HLi97KHtfezp7PHtbe317fHsxezJ7M3s0e0cuNXsgPjEuMS42ezd7OHtoJVIlUiVSJT4jOXtHLkcucn5+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IEcrKyArICAgICAgICAgIiwKIiAgICAgICsgKyA1K34gfiB+IH4gfiApI0cuRy4we34gfiBhe2J7QSxje2R7dytJLTcsNyxKLVhAZCple0cuRy5mezJAWCtnJ2cndSN9I08sZ3tHLkcuaHtMK2l7M35qe2t7bHtte257WiNve3B7Ry5xe3wufC5FO3J7c3t0e2klVCVpJT4jdXt2e0cuRy42K34gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gRysrICsgICAgICAgICAiLAoiICAgICAgKyArIDUrfiB+IH4gfiB+IHd7Ry5HLlUhfiB+IH4geHt5e3p7Sz5BezcsQns+LFAqVC40QEN7Ry5HLkR7WStYK1crViNXK1YjX0BFe0cuRy5Ge0wuTCswJXY+Mi5HezZ7SHsuPi4+cHtHLkl7SnsnLFY+aD5Le2ApTCRoJWklTHtBKi8gRy5HLk17fiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiBHKysgKyAgICAgICAgICIsCiIgICAgICArICsgNSt+IH4gfiB+IH4gfiBOe0cudS5yfn4gfiB+IDUgT3tQe14rUXtSe1N7VHtkKjRAVXtWe0cuYiFXezFAOiZ1I1crWHtZe0cuRy5HLkZ7aytqKzY7NjtqKzY7ID4xLi4+KSxwe0cuWnsrPlU+VT4vIWB7NiUyPV0kTHt+IH4gfiAgXUcuLl1+IH4gfiB+IH4gW34rXUBdI11+IH4gfiB+IH4gfiB+IEcrKyArICAgICAgICAgIiwKIiAgICAgICsgKyA1K34gfiB+IH4gfiB+ICRdfStHLlsrQF0lXSRdJl0qXT1dLV1HLkcuYiE7XVAqNEBaKz5dSC5HLkcuLF0nXSldVUAhXUcuej1HLkcuTkB+XTc7NztqKzY7e11dXV5dL10vXShdRy5fXTpdOl06XTxdW11SJTsjVSZBKn4gfiB+IHYuRy4uK31dMHtxLjwrNyt1LkcuRy4pI34gfiB+IH4gfiB+IH4gRysrICsgICAgICAgICAiLAoiICAgICAgKyArIDUrfiB+IH4gfiB+IH4gfiB8XUhARy5HLkcuRy5HLkcuRy5HLj4rMV0yXXJAdSszXXIrMkA0XSN+Ry5HLkcuRy5HLkcuNV02XUcuRy5HLlUrN100O2krOF02LUcuRy5HLkcuRy5HLkcuRy5HLkcuRy45XTBdYV1+IH4gfiB+IH4gYl1jXUcuRy5HLkcuRy5HLkcuJyYpJmRdfiB+IH4gfiB+IH4gfiBHKysgKyAgICAgICAgICIsCiIgICAgICArICsgNSt+IH4gfiB+IH4gfiB+IH4geC5lXWZdZ111Lj57MCtmXWhdaV1qXWtdUy4zXXIrWisyQDJAbF1tXSwrI341LW5db11wXXFdcl1HLkcuc110XUouID51XXZARy5HLkcuRy5HLkcuRy5HLkcuRy5HLnZdd11PI34gfiB+IH4gfiB+IHZ7eF1CJTIrLEBCJTcreV16XX4gfiB+IH4gfiB+IH4gfiB+IEcrKyArICAgICAgICAgIiwKIiAgICAgICsgKyBBXX4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gQl1DXURdSCUzO15+Xn5FXUZdR11IXUxAaiQ6QFMrUitdQChAaytrK2orSi5KLkouID4xLi4+JyxVPlM+SV0tPkA+Sl1LXUxdTV1+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gTl0rICsgICAgICAgICAiLAoiICAgICAgKyArICsgT11+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IFBdUV1SXS9+U107LlRdVV1WXXw7V11YXVldXUBdQE4rTitNK2orNjs2OyA+ID4xLiksJyxWPlY+VT4tPlpdYF0gXi5eK15+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiBAXisgKyAgICAgICAgICIsCiIgICAgICArICsgKyAjXn4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+ICReJV4mXipePV4tXjtePl4sXideKV4hXjZ7fl5LLksuaStLLjU7ISwuPiksJyxVPkldcj57Xl1eXl4vXihefiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiBfXjpeKyArICAgICAgICAgIiwKIiAgICAgICAgKyArIDxefiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IGZAW159XnxeMV4yXjNeNF41XjZeN14oLVojZz4gPkouID4uPi4+JyxVPjw7OF45XjBeRiRhXnIgfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IGJeKyArICsgICAgICAgICAiLAoiICAgICAgICArICsgY15kXn4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IG4gZV5mXmdeaF5pXmpea15sXm1ebl5vXnBeMTtxXnJec150XnVedl53XnheciBuIH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4geV4rICsgICAgICAgICAgICIsCiIgICAgICAgICsgKyArIHpefiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IEFeQl5DXkReRV5GXkdeSF5JXkpeS15MXihefiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gTV5OXisgKyAgICAgICAgICAgIiwKIiAgICAgICAgICArICsgT15QXn4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gNCA0IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiBHKysgKyAuICAgICAgICAgICAiLAoiICAgICAgICAgICcgKyArIFFeQi5+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiBSXisgKyArICAgICAgICAgICAgICIsCiIgICAgICAgICAgICArICsgKyBTXjIgfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiBUXlVeKyArICAgICAgICAgICAgICAgIiwKIiAgICAgICAgICAgICAgKyArICsgVl5XXn4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gWF5RICsgKyArICsgICAgICAgICAgICAgICAiLAoiICAgICAgICAgICAgICAgICsgKyArIEMuUiBZXn4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IH4gfiB+IDIgel5PXisgKyArICsgICAgICAgICAgICAgICAgICIsCiIgICAgICAgICAgICAgICAgICArICsgKyArIFsgWl5gXiAvLi8rLysvKy8rLysvKy8rLysvKy8rLysvKy8rLysvKy8rLysvKy8rLysvKy8rLysvKy8rLysvKy8rLysvKy8rLysvKy8rLysvKy8rLysvKy8rLysvKy8rLysvKy8rLysvKy8rLysvKy8rLysvKy8rLysvKy8rLysvQC8jLyQvJS8rICsgKyArICAgICAgICAgICAgICAgICAgICAgIiwKIiAgICAgICAgICAgICAgICAgICAgLiArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICAgICAgICAgICAgICAgICAgICAgICAiLAoiICAgICAgICAgICAgICAgICAgICAgICAgLiArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyAgICAgICAgICAgICAgICAgICAgICAgICAgICIsCiIgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyArICsgKyAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgIiwKIiAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAiLAoiICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICAgICJ9Owo="

#################################################################
#
# Main subroutines
#
#################################################################

# Prepare to plot...
app = PyQt5.QtWidgets.QApplication(sys.argv)	
form = plotEsg()
app.exec_()
