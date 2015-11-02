#!/usr/bin/env python
"""
 wrapper script for running Elprep application
"""

import sys, argparse, os, tempfile, subprocess, shutil

# Cleanup tmp_dir
def cleanup_before_exit( tmp_dir ):
   if tmp_dir and os.path.exists( tmp_dir ):
        shutil.rmtree( tmp_dir )


# Function to create a symbolic link of the input.dat file to a /tmpDir/input.sam or /tmpDir/input.bam file
def symlink_elprep_input_filename( filename, filetype, target_dir = None ):
    if target_dir is None:
        target_dir = os.getcwd()
    print filename
    base_filename = os.path.basename(filename)
    print base_filename
    new_filename = os.path.join( "%s.%s" % (base_filename, filetype ) )
    new_input_file = os.path.join ( "%s/%s" % (target_dir, new_filename) )
    print new_input_file
    #Adding a check for file existent or not
    if not os.path.isfile(new_input_file):
        shutil.copyfile ( filename, new_input_file)
        #shutil.copy2( filename, new_input_file)
	#os.symlink( filename, new_input_file )
    return new_input_file


def __main__():
    #Parse Command Line
    parser = argparse.ArgumentParser()

    parser.add_argument( '-i', '--input', dest='input', action='store', type=str, help='Input SAM/BAM file' )
    parser.add_argument( '-f', '--input_file_type', dest='input_file_type', action='store', help='Input file type for SAM/BAM file' )
    parser.add_argument( '-o', '--output', dest='output', help='The output BAM dataset' )
    parser.add_argument( '-p', '--pass_through', dest='pass_through_options', action='store', type=str, help='These options are passed through directly to the application without any modification.' )

    options = parser.parse_args()
    
    # Create a tmpdir for storing input and output files with the file extensions that Elprep prefers
    tmp_dir = tempfile.mkdtemp( prefix='tmp-elprep-' )
    print tmp_dir

    if options.input and options.input_file_type:
         tmp_input_file = tempfile.NamedTemporaryFile( dir=tmp_dir )
         tmp_input_filename = tmp_input_file.name
         tmp_input_file.close()
         new_tmp_input = os.path.join( "%s.%s" % (tmp_input_filename, options.input_file_type ) )
         shutil.copyfile ( options.input, new_tmp_input)
         #new_tmp_input = symlink_elprep_input_filename(options.input, options.input_file_type, tmp_input_filename, target_dir = tmp_dir)
	 print new_tmp_input 

    if options.output:
         tmp_output_file = tempfile.NamedTemporaryFile( dir=tmp_dir )
         tmp_output_filename = tmp_output_file.name
         tmp_output_file.close()
         tmp_output_filename = os.path.join( "%s.bam" % (tmp_output_filename ) )
         print tmp_output_filename

    #Construct the new command line 
    cmd = 'elprep %s %s %s' % ( options.pass_through_options, new_tmp_input, tmp_output_filename )
    #cmd = 'elprep %s %s %s' % ( options.pass_through_options, tmp_input_filename, tmp_output_filename )
    #cmd = 'elprep %s /home/aprabh2/tmp.bam %s' % ( tmp_input_filename, options.pass_through_options )
    print "Command:"
    print cmd
 
    #Run the Elprep program
    try:
        subprocess.check_output(args=cmd, stderr=subprocess.STDOUT, shell=True)
        shutil.move( tmp_output_filename, options.output )
    except subprocess.CalledProcessError, e:
        print "!!!!!!!!!!!! Elprep ERROR: stdout output:\n", e.output

    cleanup_before_exit( tmp_dir )

if __name__=="__main__": __main__()
