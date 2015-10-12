#!/usr/bin/env python
# -*- coding: utf-8 -*- 

from PIL import Image
from PIL import ImageEnhance
from sys import argv
import os
import optparse

parser = optparse.OptionParser()
parser.add_option('-i', action="append", dest="image")
parser.add_option('-o', action="store", dest="output")
parser.add_option('-b', action="store", dest="id")
parser.add_option('-d', action="store", dest="directory")
parser.add_option('-e', action="store", dest="enhancement")
parser.add_option('-j', action="store", dest="value")
(options, args) = parser.parse_args()

def tool_Tester(image, output, id, directory, enhancement, value):
	for i,item in enumerate(options.image):
		if enhancement == "Color":
			if i == 0:
				myimage = Image.open(item)
				enhancer = ImageEnhance.Color(myimage)
				myimage = enhancer.enhance(float(value))	
				output_image = myimage.save(output, format="png")
		
			else:
				myimage = Image.open(item)
				enhancer = ImageEnhance.Color(myimage)
				myimage = enhancer.enhance(float(value))
				output_image =  myimage.save("primary_"+id+"_output" + str(i+1) + "_visible_png", format="png")
		if enhancement == "Brightness":
			if i == 0:
                                myimage = Image.open(item)
                                enhancer = ImageEnhance.Brightness(myimage)
                                myimage = enhancer.enhance(float(value))
                                output_image = myimage.save(output, format="png")

                        else:
                                myimage = Image.open(item)
                                enhancer = ImageEnhance.Brightness(myimage)
                                myimage = enhancer.enhance(float(value))
                                output_image =  myimage.save("primary_"+id+"_output" + str(i+1) + "_visible_png", format="png")

		if enhancement == "Contrast":
			if i == 0:
                                myimage = Image.open(item)
                                enhancer = ImageEnhance.Contrast(myimage)
                                myimage = enhancer.enhance(float(value))
                                output_image = myimage.save(output, format="png")

                        else:
                                myimage = Image.open(item)
                                enhancer = ImageEnhance.Contrast(myimage)
                                myimage = enhancer.enhance(float(value))
                                output_image =  myimage.save("primary_"+id+"_output" + str(i+1) + "_visible_png", format="png")


		if enhancement == "Sharpness":
			if i == 0:
                                myimage = Image.open(item)
                                enhancer = ImageEnhance.Sharpness(myimage)
                                myimage = enhancer.enhance(float(value))
                                output_image = myimage.save(output, format="png")

                        else:
                                myimage = Image.open(item)
                                enhancer = ImageEnhance.Sharpness(myimage)
                                myimage = enhancer.enhance(float(value))
                                output_image =  myimage.save("primary_"+id+"_output" + str(i+1) + "_visible_png", format="png")

	return output_image
	
tool_Tester(options.image, options.output, options.id, options.directory, options.enhancement, options.value)
