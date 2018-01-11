import sys
import os
import Gnuplot
import Numeric
import math
from numpy import *
from numpy.fft import *
from numpy.linalg import *
from numpy.oldnumeric.linear_algebra import *
from Numeric import *
from CubicSpline import *

""" 
 This script reads specified directory and creates a Bubble Report.: 

   1. Converts the points to there polar coordinates. 
   2. Fits a cubic spline to the converted vectors. 
   3. Finds an equal spaced set of points around the cubic spline found.
   4. Takes the FFT of equaly spaced points (equal spacing in theta). 
   5. Plots a specified mode coefficient of the FFT vs. Time. 

   To Run:
   [user@machine directory]$ python FFTBubbleIBpoint.py 

"""

# ============================================================================ #
# Run Time options                                                             #
# ============================================================================ #
DataType         = "vtk"              # Options are "vtk" or "dat"
integerData      = 1                  # integer style data sets or 0 for old time series style

showonscreen     = 0                  # Set to 1 to show plots as made to screen 
saveFFTPlots     = 0
saveBubblePlots  = 0
NoPlots          = 0

Directory        = "./RunOn3_23_2010Bree/"
DirectoryOut     = Directory
DirectoryShort   = Directory

rate = 0.95
frq  = 4.45

FileBase         = "BubblePts_time_"

matlabfile       = "FFTData.m"        # Name of the matlab file to write the coefficents to. 
latexfileName    = "BubbleReport.tex" # Name of the matlab file to write the coefficients to. 
FirstAndLastName = "NotGiven";
FFTPlotName      = "NotGiven";
SecondModeName   = "NotGiven";

sys.path.append(Directory)
from PostProcessParams import *

# Time:                    
#T         = 1.25
#dt        = 0.0005 
#timeSteps = 2500           
#numout    = 25 

# ============================================================================ #
# Local Variables                                                              #
# ============================================================================ #
centerx = (xmax + xmin)/2.0; centery = (ymax + ymin)/2.0; 
IB_r_th_File_name = 'gnuplot_IB_r_th_file.txt'
IB_r_th_File  = open(IB_r_th_File_name,"w")
pi            = 3.14159265358979
Time          = 0.0
dt            = T/timeSteps
Arc_Length    = 0.0                      #  Total Arclength of splined bubble (use a lot of points)
npts_spline   = 512                      #  Number of points on the spline (power of 2 for FFT) 
a_lim         = 0.0                      #  Starting point of arch
b_lim         = 2.0*pi                   #  Ending point of arch
maxItr        = 600                      # Maximum number of Newton Iterations to take per data set. 
errtol        = 1e-14                    # Stop the Newton Iterations if the diff in solution is less than this value.
FFTScaling    = (2.0/(npts_spline))*(1.0/(bub_A))
dth           = (b_lim - a_lim)/npts_spline
pointCounter  = 0    
saved_index   = 0      
latexfile     = DirectoryOut + latexfileName
LATEX         = open(latexfile,"w")
totalFrames   = timeSteps/numout + 1
frameStep     = T/(totalFrames-1)
dt            = T/timeSteps
pointCounter  = 0
FileIndex     = 0
ReadingPoint  = 0

# SOME ARRAYS
IBpoints_x    = zeros(totalFrames*bub_nIBpts, Float); IBpoints_x.shape=(bub_nIBpts,totalFrames)
IBpoints_y    = zeros(totalFrames*bub_nIBpts, Float); IBpoints_y.shape=(bub_nIBpts,totalFrames)
IBpoints_r    = zeros(totalFrames*bub_nIBpts, Float); IBpoints_r.shape=(bub_nIBpts,totalFrames)
IBpoints_th   = zeros(totalFrames*bub_nIBpts, Float); IBpoints_th.shape=(bub_nIBpts,totalFrames)
time          = zeros(totalFrames, Float);           time.shape=(totalFrames,1)
spline_theta  = zeros(npts_spline, Float);   spline_theta.shape=(npts_spline,1)
spline_radius = zeros(npts_spline, Float);   spline_radius.shape=(npts_spline,1)
fft_vector    = zeros(npts_spline, Float);   fft_vector.shape=(npts_spline,1) 
mode          = zeros(npts_spline,Float);    mode.shape=(npts_spline,1)
amp_pmode     = zeros(totalFrames,Float);        amp_pmode.shape=(totalFrames,1)
amp_2pmode    = zeros(totalFrames,Float);        amp_2pmode.shape=(totalFrames,1)
ana_pmode     = zeros(totalFrames,Float);        ana_pmode.shape=(totalFrames,1)
newTime       = zeros(totalFrames,Float);        newTime.shape=(totalFrames,1)

frameStep   = T/totalFrames

yRange = 'set yrange [' + str(-1.5*bub_R) + ':' + str(1.5*bub_R) + ']'
xRange = 'set xrange [' + str(-1.5*bub_R) + ':' + str(1.5*bub_R) + ']'

if(DataType == "dat"):
    for filecount in range(totalFrames):   
        TIME = frameStep*filecount
        if(integerData == 1):
            datFile = Directory + FileBase + str('%(#)08d.dat' % {"#":filecount})
        if(integerData == 0):
            datFile = Directory + FileBase + str('%(#)0.6f.dat' % {"#":TIME})
        myfile  = open(datFile,"r")
        print "======================================================"
        print "Reading data from " + datFile 
        print "======================================================"
        #read line into array
        for line in myfile.readlines():
            numbers = map(float, line.split())
            IBpoints_x[pointCounter,saved_index] = float(numbers[0])
            IBpoints_y[pointCounter,saved_index] = float(numbers[1]) 
            pointCounter = pointCounter + 1
            if(pointCounter == bub_nIBpts):
                pointCounter = 0
                saved_index = saved_index + 1

print "-------------------------totalFrames " , totalFrames
if(DataType == "vtk"):
    for filecount in range(totalFrames):   
        TIME = frameStep*filecount
        datFile = Directory + FileBase + str('%(#)08d.vtk' % {"#":filecount})
        myfile  = open(datFile,"r")
        #print "==============================================================================================================================================="
        #print "Reading data from " + datFile 
        #print "===============================================================================================================================================\n"
        #read line into array
        for line in myfile.readlines():
            numbers = map(str, line.split())
            if((ReadingPoint == 1) and (pointCounter < bub_nIBpts)):
                #print "reading in point ", pointCounter, "x = ", numbers[0], " y = ", numbers[1]
                IBpoints_x[pointCounter,FileIndex] = float(numbers[0])
                IBpoints_y[pointCounter,FileIndex] = float(numbers[1])
                pointCounter = pointCounter + 1
                if(pointCounter == nIBpts):
                    ReadingPoint = 0
            if(ReadingPoint == 0):
                if(len(numbers) > 0):
                    if(numbers[0] == "POINTS"):
                        nIBpts = int(numbers[1])
                        #print "The number of IB points is", nIBpts
                        ReadingPoint = 1
                        pointCounter = 0
            if(pointCounter == bub_nIBpts):
                pointCounter = 0
                ReadingPoint = 0
                FileIndex = FileIndex + 1

totalSaves = totalFrames
i = 0
while(i < totalSaves):
   time[i,0]=(i*numout*dt)
   i=i+1

#============================================================
print "Calculating the IB points in terms of radius and theta" 
print "======================================================"
i = 0
j = 0
for j in range(totalSaves):
   for i in range(bub_nIBpts):
      IBpoints_r[i,j] = ((IBpoints_x[i,j]-centerx)**2.0 + (IBpoints_y[i,j]-centery)**2.0)**(0.5)
      if((IBpoints_y[i,j]-centery) < 0):
         IBpoints_th[i,j] = 2.0*pi - arccos((IBpoints_x[i,j]-centerx)/IBpoints_r[i,j])
      else: 
         IBpoints_th[i,j] = arccos((IBpoints_x[i,j]-centerx)/IBpoints_r[i,j])

      IB_r_th_File.write(str(IBpoints_r[i,j]) + "  " + str(IBpoints_th[i,j]) + "\n") 
   IB_r_th_File.write("\n\n") 

# Set up an equally spaced theta vector. 
for i in range(len(spline_theta)):
   spline_theta[i,0] = float((float(i))/float(len(spline_theta)))*2.0*pi
 
for k in range(totalSaves):
   #============================================================
   print "Calculating a Cubic spline for the data set " + str(k) 
   print "======================================================"

   #Spline = cubic_spline(xvec=IBpoints_th[:,k],yvec=IBpoints_r[:,k])
   Spline = fast_cubic_spline(xvec=IBpoints_th[:,k],yvec=IBpoints_r[:,k])
   coefficents = Spline.periodic()
   
   #============================================================
   #print "Calculate set of Equaly spaced points around bubble   " 
   #print "======================================================"
   # FIRST CALCULATE THE TOTAL ARCLENGTH OF THE BUBBLE ON A GIVEN 
   # INTERVAL USING THE CUBIC SPLINE FOUND ABOVE. 

   j = 0 
   for i in range(npts_spline):
      # Calculate the arch lenght of the spline using points equally spaced in theta. 
      foundInterval = 'false'
      while(foundInterval == 'false'):
         
         ################## USE THIS BLOCK WITH PLAIN "cubic_spline" #########################
         #if(IBpoints_th[(j-1)%bub_nIBpts,k] <= spline_theta[i,0] < IBpoints_th[j%bub_nIBpts,k]):
         #
         #   spline_radius[i,0] =  (coefficents[(j-1)% bub_nIBpts,0]
         #                        + coefficents[(j-1)% bub_nIBpts,1]*( spline_theta[i,0]) 
         #                        + coefficents[(j-1)% bub_nIBpts,2]*((spline_theta[i,0])**2.0)  
         #                        + coefficents[(j-1)% bub_nIBpts,3]*((spline_theta[i,0])**3.0))
         #   foundInterval = 'true'    
         #   mode[i] = i
         #
         #if( (IBpoints_th[(j-1)%bub_nIBpts,k] > IBpoints_th[j%bub_nIBpts,k]) and (0.0 <= spline_theta[i,0] < IBpoints_th[j%bub_nIBpts,k])):
         #   spline_radius[i,0] =  (coefficents[(j-1)% bub_nIBpts,0]
         #                        + coefficents[(j-1)% bub_nIBpts,1]*( spline_theta[i,0] + 2.0*pi) 
         #                        + coefficents[(j-1)% bub_nIBpts,2]*((spline_theta[i,0] + 2.0*pi)**2.0)  
         #                        + coefficents[(j-1)% bub_nIBpts,3]*((spline_theta[i,0] + 2.0*pi)**3.0))
         #   foundInterval = 'true'    
         #   mode[i] = i
         #
         #if((IBpoints_th[(j-1)%bub_nIBpts,k] > IBpoints_th[j%bub_nIBpts,k]) and (IBpoints_th[(j-1)%bub_nIBpts,k] <= spline_theta[i,0])):
         #   spline_radius[i,0] =  (coefficents[(j-1)% bub_nIBpts,0]
         #                        + coefficents[(j-1)% bub_nIBpts,1]*( spline_theta[i,0]) 
         #                        + coefficents[(j-1)% bub_nIBpts,2]*((spline_theta[i,0])**2.0)  
         #                        + coefficents[(j-1)% bub_nIBpts,3]*((spline_theta[i,0])**3.0))
         #   foundInterval = 'true'    
         #   mode[i] = i
         #else:
         #   j = j+1
         ##################################################################################END

         ################## USE THIS BLOCK WITH "fast_cubic_spline" ##########################
         if(IBpoints_th[(j-1)%bub_nIBpts,k] <= spline_theta[i,0] < IBpoints_th[j%bub_nIBpts,k]):

            spline_radius[i,0] =  (coefficents[(j-1)% bub_nIBpts,0]
                                 + coefficents[(j-1)% bub_nIBpts,1]*( spline_theta[i,0]-IBpoints_th[(j-1)%bub_nIBpts,k]) 
                                 + coefficents[(j-1)% bub_nIBpts,2]*((spline_theta[i,0]-IBpoints_th[(j-1)%bub_nIBpts,k])**2.0)  
                                 + coefficents[(j-1)% bub_nIBpts,3]*((spline_theta[i,0]-IBpoints_th[(j-1)%bub_nIBpts,k])**3.0))
            foundInterval = 'true'    
            mode[i] = i

         if( (IBpoints_th[(j-1)%bub_nIBpts,k] > IBpoints_th[j%bub_nIBpts,k]) and (0.0 <= spline_theta[i,0] < IBpoints_th[j%bub_nIBpts,k])):
            spline_radius[i,0] =  (coefficents[(j-1)% bub_nIBpts,0]
                                 + coefficents[(j-1)% bub_nIBpts,1]*( spline_theta[i,0]-IBpoints_th[(j-1)%bub_nIBpts,k] + 2.0*pi) 
                                 + coefficents[(j-1)% bub_nIBpts,2]*((spline_theta[i,0]-IBpoints_th[(j-1)%bub_nIBpts,k] + 2.0*pi)**2.0)  
                                 + coefficents[(j-1)% bub_nIBpts,3]*((spline_theta[i,0]-IBpoints_th[(j-1)%bub_nIBpts,k] + 2.0*pi)**3.0))
            foundInterval = 'true'    
            mode[i] = i

         if((IBpoints_th[(j-1)%bub_nIBpts,k] > IBpoints_th[j%bub_nIBpts,k]) and (IBpoints_th[(j-1)%bub_nIBpts,k] <= spline_theta[i,0])):
            spline_radius[i,0] =  (coefficents[(j-1)% bub_nIBpts,0]
                                 + coefficents[(j-1)% bub_nIBpts,1]*( spline_theta[i,0]-IBpoints_th[(j-1)%bub_nIBpts,k]) 
                                 + coefficents[(j-1)% bub_nIBpts,2]*((spline_theta[i,0]-IBpoints_th[(j-1)%bub_nIBpts,k])**2.0)  
                                 + coefficents[(j-1)% bub_nIBpts,3]*((spline_theta[i,0]-IBpoints_th[(j-1)%bub_nIBpts,k])**3.0))
            foundInterval = 'true'    
            mode[i] = i
         else:
            j = j+1
         ##################################################################################END

   for i in range(len(spline_radius)):
     fft_vector[i,0] = FFTScaling*(spline_radius[i,0] - bub_R)

   # Walk around the spline and populate an equally spaced set of points. 
   f_fft   = fft(fft_vector[:,0])
   f_fft_r = fft(fft_vector[:,0]).real
   f_fft_i = fft(fft_vector[:,0]).imag
   amp_pmode[k,0] = f_fft_r[int(bub_p)]
   amp_2pmode[k,0] = f_fft_r[int(2*bub_p)]

   # ========================================================================= #
   # Plot the functions fft (real part) using gnuplot.                         #
   # ========================================================================= #

   if(saveFFTPlots == 1):
       if(k%2 == 0):
           SavePlotName = DirectoryOut + str('Function_FFT_real_t') + str('%(#)0.6f.eps' % {"#":time[k,0]}) 
           g = Gnuplot.Gnuplot(debug=1)
           d = Gnuplot.Data(mode[:,0],f_fft_r)
           Title = 'Plot of Real Part of FFT(r(theta)) Time = '+ str(time[k,0])
           xaxis = 'mode'
           yaxis = 'FFT(r(theta))' 
           g.title(Title)     
           g('set data style lines') # give gnuplot an arbitrary command
           yFFt = 'set yrange [0.2:' + str(1.2*bub_R) + ']'
           g(yFFt)
           #g('set xrange [0:1.01]')
           g.xlabel(xaxis)
           g.ylabel(yaxis)
           if(showonscreen == 1):
               g.plot(d)
               raw_input('Please press return to continue...\n')

           if(SavePlotName != ''):
               outFile = 'set output \"' + SavePlotName + '\"'
               g(outFile)
               g('set term postscript eps enhanced color solid')
               g.plot(d)

   # ========================================================================= #
   # Plot the Spline and the immersed boundary points at current time          #
   # ========================================================================= # 
   if(saveBubblePlots == 1):
       SavePlotName = DirectoryOut + str('SplineandIBpoints_t') + str('%(#)0.6f.eps' % {"#":time[k,0]}) 
       g = Gnuplot.Gnuplot(debug=1)
       d = Gnuplot.Data(IBpoints_th[:,k],IBpoints_r[:,k])
       d2 = Gnuplot.Data(spline_theta[:,0],spline_radius[:,0])
       Title = 'Bubble'
       xaxis = 'Plot of IB points and Spline Time = ' + str(time[k,0])
       yaxis = '' 
       g.title(Title)     
       g('set data style lines') # give gnuplot an arbitrary command
       g('set polar')
       g('set size square 1,1')       
       g(xRange)
       g(yRange)
       g.xlabel(xaxis)
       g.ylabel(yaxis)
       if(showonscreen == 1):
           g.plot(d,d2)
           raw_input('Please press return to continue...\n')

       if(SavePlotName != ''):
           outFile = 'set output \"' + SavePlotName + '\"'
           g(outFile)
           g('set term postscript eps enhanced color solid')
           g.plot(d,d2)


# ========================================================================
# Fit the function y = exp(-rate*t)*cos(frequency*t) to the fftmode         
# ========================================================================
npts    = totalSaves;
J       = ones(npts*2, Float); J.shape=(npts,2)
JTJ     = ones(4, Float);      JTJ.shape=(2,2)
ratefrq = ones(2, Float);      ratefrq.shape=(2,1)
r       = ones(npts, Float);   r.shape=(npts,1)
JTr     = zeros(2, Float);     JTr.shape=(2,1)
# Give an initial guess at the decay rate and the frequency. 
# (rate = 1.5 frq  = 13.0 : Re = 10, R = 1, A = 0.2, sigma_b = 25)

norm_xiandipo = 1.0
iteration = 0

while((norm_xiandipo > errtol) and (iteration < maxItr)):
    iteration = iteration + 1
    # Populate the Jacobian
    for i in range(npts):
        x_i = time[i,0]
        #print x_i
        J[i,0] = x_i*cos(frq*x_i)*exp(-rate*x_i) # dy/d(rate)
        J[i,1] = x_i*sin(frq*x_i)*exp(-rate*x_i) # dy/d(frq)
        r[i,0] = amp_pmode[i,0] - exp(-rate*x_i)*cos(frq*x_i)
    
    # Set up a least squares system 
    JTJ = matrixmultiply(transpose(J),J)
    rtmp = transpose(matrixmultiply(transpose(J),r[:,0]))

    JTr[0,0] = -rtmp[0]
    JTr[1,0] = -rtmp[1]

    System = gmres(Matrix=JTJ,RHS=JTr,x=ratefrq,Tol=1e-14,maxIts=len(JTr))
    coefficents, error, totalIters = System.solve()

    # Shift Solution values
    rateold = rate; frqold = frq;           
    rate = rate + coefficents[0,0];     
    frq = frq + coefficents[1,0];     
  
    ratefrq[0,0] = rate; 
    ratefrq[1,0] = frq;  

    # calculate the norm. 
    string = "iteration " + str(iteration) + " rate = " + str(rate) + " frq = " + str(frq) +"\n"
    print string
    # REPORT.write(string)  
    norm_xiandipo = ((rateold-rate)**2.0 + (frqold-frq)**2.0)**0.5

    # end Solve of least squares problem.

#print "The decay rate and frequence found by iteration ", iteration, " are rate = ", rate, " frequency = ", frq, "\n" 

if(NoPlots != 1):
    # ============================================================================ #
    # Plot the Amplitude of the p mode in the fft of the given data set            #
    # ============================================================================ #

    # Set up a comparison function
    for i in range(size(amp_pmode[:,0])):
        t = 1.0*time[i,0]
        t2 = time[i,0]
        newTime[i,0] = t
        frequency = ((bub_p**3 - bub_p)*(bub_sprKlink/(bub_R**3.0)))**(0.5)
        if(bub_p == 3.0):
            #ana_pmode[i,0] = bub_R*cos(t2*((24.0*pi)**(0.5)))
            #ana_pmode[i,0] = bub_R*cos(t2*frequency) 
            ana_pmode[i,0] = exp(-rate*t)*cos(frq*t)
        if(bub_p == 2.0):
            ana_pmode[i,0] = bub_R*cos(t2*((6.0*pi)**(0.5)))
    FFTPlotName = DirectoryOut + str('AmpModePlot') + str(int(bub_p))
    SavePlotName =  FFTPlotName +'.eps'
    g  = Gnuplot.Gnuplot(debug=1)
    d  = Gnuplot.Data(time[:,0],amp_pmode[:,0])
    d2 = Gnuplot.Data(newTime[:,0],ana_pmode[:,0])
    Title = 'Amplitude of p=' + str(int(bub_p)) + ' FFT mode'
    xaxis = 'time'
    yaxis = 'Amplitude' 
    g.title(Title)     
    g('set data style lines') # give gnuplot an arbitrary command
    #g('set yrange [0:1.0]')
    #g('set xrange [0:1.01]')
    g.xlabel(xaxis)
    g.ylabel(yaxis)
    if(showonscreen == 1):
        g.plot(d,d2)
        raw_input('Please press return to continue...\n')

    if(SavePlotName != ''):
        outFile = 'set output \"' + SavePlotName + '\"'
        g(outFile);
        string = 'set terminal postscript eps enhanced color solid';
        g(string);
        g.plot(d,d2)

    # ============================================================================ #
    # Plot the Amplitude of the 2p mode in the fft of the given data set            #
    # ============================================================================ #

    # Set up a comparison function
    for i in range(size(amp_pmode[:,0])): 
        t = 1.0*time[i,0]
        t2 = time[i,0]
        newTime[i,0] = t
        if(bub_p == 3.0):
            ana_pmode[i,0] = bub_A*((43.0/280.0) + (59.0/152.0)*cos(t2*2.0*((24.0*pi)**(0.5))) - (1441.0/2660.0)*cos(((210.0*pi)**(0.5))*t2) )
        if(bub_p == 2.0): 
            ana_pmode[i,0] = bub_A*((3.0/20.0) + (1.0/3.0)*cos(t2*2.0*((6.0*pi)**(0.5))) - (29.0/60.0)*cos(((60.0*pi)**(0.5))*t2) )

    SecondModeName = DirectoryOut + str('Amp2ModePlot') + str(int(bub_p)) 
    SavePlotName   = SecondModeName + '.eps'
    g  = Gnuplot.Gnuplot(debug=1)
    d  = Gnuplot.Data(time[:,0],amp_2pmode[:,0])
    d2 = Gnuplot.Data(newTime[:,0],ana_pmode[:,0])

    Title = 'Amplitude of p=' + str(int(2*bub_p)) + ' FFT mode'
    xaxis = 'time'
    yaxis = 'Amplitude' 
    g.title(Title)     
    g('set data style lines') # give gnuplot an arbitrary command
    #g('set yrange [0:1.0]')
    #g('set xrange [0:1.01]')
    g.xlabel(xaxis)
    g.ylabel(yaxis)
    if(showonscreen == 1):
        g.plot(d)
        raw_input('Please press return to continue...\n')

    if(SavePlotName != ''):
        outFile = 'set output \"' + SavePlotName + '\"'
        g(outFile)
        g('set terminal postscript eps enhanced color solid')
        g.plot(d)

    # ============================================================================ #
    # Plot of the first and last IB bubble using gnuplot.                          #
    # ============================================================================ #
    FirstAndLastName = str('FirstAndLastIBPoints') + str(int(bub_p)) 
    SavePlotName = DirectoryOut + FirstAndLastName + str('.eps')
    g = Gnuplot.Gnuplot(debug=1)
    d = Gnuplot.Data(IBpoints_th[:,0],IBpoints_r[:,0])
    d2 = Gnuplot.Data(IBpoints_th[:,(totalFrames-1)],IBpoints_r[:,(totalFrames-1)])
    Title = 'Plot of Bubble IB points at Time = ' + str(time[0,0]) + ' in red, and T = ' + str(time[(totalFrames-1),0]) + " in green"
    xaxis = 'Bubbles'
    yaxis = '' 
    g.title(Title)     
    g('set data style lines') # give gnuplot an arbitrary command
    g('set polar')
    #g('set size square 1,1')
    g(yRange)
    g(xRange)
    g.xlabel(xaxis)
    g.ylabel(yaxis)
    if(showonscreen == 1):
        g.plot(d,d2)            
        g('set terminal x11')
        raw_input('Please press return to continue...\n')

    if(SavePlotName != ''):
        outFile = 'set output \"' + SavePlotName + '\"'
        g(outFile)
        g('set terminal postscript eps enhanced color solid')
        g.plot(d,d2)

# =========================================================================
# Write LaTeX file for further use of the velocity values calculated. 
# =========================================================================
LATEX.write("\\documentclass[11pt]{article} \n")
LATEX.write("\\usepackage{geometry} \\geometry{hmargin={1in,1in},vmargin={1in,1in}} \n")
LATEX.write("\\usepackage{amsmath,amsfonts,amssymb,color} \n")
LATEX.write("\\usepackage[dvips]{graphics} \n")
LATEX.write("\\usepackage[dvips]{graphicx} \n")
LATEX.write("\\begin{document} \n")
LATEX.write("\\begin{center}\n")
LATEX.write("{\\bf \\Large Bubble Report} \n")
LATEX.write("\\end{center}\n")
# Section Parameters
LATEX.write("\\section{Parameters:} \n")
LATEX.write("This is a report that examines an immersed fiber initially expressed in terms of cosine perterbations of a circle. The initial conditions are set using:\n")
string = "\\begin{equation} \n r = R + A \\cos( n\\theta) \\label{bubble}\n\\end{equation} \n"; LATEX.write(string);
string = "\n\n \\vskip.2in \\noindent "; LATEX.write(string); 
string = "This code was run on " + str(Date) + " with the parameters given in the following table:\n "; LATEX.write(string);
string = "\\begin{table}[!htp] \n"; LATEX.write(string)
string = "\\begin{center} \n"; LATEX.write(string)
string = "\\vskip.1in \n"; LATEX.write(string)
string = "\\begin{tabular}{|l|l|l|l|l|} \\hline \n"; LATEX.write(string)
string = "Grid: &  Time: & Fluid: & IB: & Bubble: \\\\  \\hline \\hline \n"; LATEX.write(string);
string = "$[x_{min} , x_{max}] $ =  ["+str(xmin)+","+str(xmax)+"] & $T$= "+str(T)+" & $Re$ = "+str(Re)+"& npts = "+str(nIBpts)+"& $R $ = "+ str(bub_R) + "\\\\ \n"; LATEX.write(string); 
string = "$[y_{min} , y_{max}] $ =  ["+str(ymin)+","+str(ymax)+"] & $N$= "+str(timeSteps) + "& $\\lambda$ = " + str(We) + " & $\\sigma$ = " + str(bub_sprKlink) + "& $A$ = " 
string = string + str(bub_A) + "\\\\ \n"; LATEX.write(string); 
string = "$(x_{pts},y_{pts})$    =  ("+str(xpts)+","+str(ypts)+") & $\\mbox{total}_{\mbox{saves}}$ = "
string = string + str(timeSteps/numout)+" & $\\alpha$ = " + str(alpha) + " & & $n$ = " + str(bub_p)  + "\\\\ \n"; LATEX.write(string);         
string = "& $\\Delta t$ = " + str(T/timeSteps) + "& $\\epsilon$ = " + str(epsilon) +" & & \\\\ \n"; LATEX.write(string); 
string = "&  & $\\nu$ = " + str(nu) + " &  & \\\\ \\hline \n "; LATEX.write(string); 
string = "\\end{tabular} \n"; LATEX.write(string); 
string = "\\caption{Numerical Model Parameters} \n"; LATEX.write(string)
string = "\\end{center}\n"; LATEX.write(string);
string = "\\label{parameters} \n"; LATEX.write(string); 
string = "\\end{table} \n\n \\noindent "; LATEX.write(string); 
string = "The Giesekus constitutive model for viscoelastic fluid flow was used: \n" 
string = string + "\\newcommand{\\bu}{\\mbox{\\boldmath $u$}}"
string = string + "\\newcommand{\\bs}{\\mbox{\\boldmath ${\\sigma}$}}"
string = string + "\\begin{eqnarray*}\n "
string = string + "\\bs + \\lambda \\left( {\\color{black}\\frac{\\partial  \\bs}{\\partial t }}+ \\bu \\cdot \\nabla \\bs -(\\nabla\\bu) \\bs - \\bs(\\nabla\\bu)^T " 
string = string + "- \\nu\Delta \\bs + \\epsilon \\bs^2  \\right) - 2\\alpha {\\bf d}(\\bu) & = & 0 \\;\\; \\mbox{in } \\Omega \\\\ \n" 
string = string + "Re\\left({\\color{black}\\frac{\\partial \\bu}{\\partial t}} + \\bu\\cdot\\nabla \\bu\\right) + \\nabla p "
string = string + "- 2(1-\\alpha)\\nabla \\cdot {\\bf d}(\\bu) - \\nabla \\cdot \\bs & = & {\\bf f} + {\\bf F_{IB}} \;\; \mbox{in } \Omega \\\\ \n"
string = string + "\\nabla \\cdot \\bu & = & 0  \\;\\; \\text{in} \\;\\; \\Omega \n"  
string = string + "\\end{eqnarray*} \n\n"; LATEX.write(string);
# Forces on the boundary:
string = "\\section{Immersed boundary forces:}\n"
string = string + "The Forces are set based on the curvature of the immersed fiber. With the curve defined in polar coordinates "
string = string + " as $ r(\\theta)$ the forces at point $\\theta$ on the fiber are: \n"
string = string + "\\[ F(\\theta) = \\sigma \\kappa(\\theta) {\\bf n}  \\]\n where \n "
string = string + "\\[ \\kappa(\\theta) = \\frac{(r^2 + 2(r')^2 - r r'')}{\\left((r)^2 + (r')^2 \\right)^{\\frac{3}{2}}} \\]\n "
string = string + "is the curvature of the fiber, ${\\bf n}$ is the outward normal and $\\sigma$ is a material parameter. "  
string = string + "The values of $r, r'$ and $r''$ are the radius, and its derivates with respect to $\\theta$."; LATEX.write(string); 
# Plot of the immersed boundary
string = "\\section{Plot of immersed boundary: }\n\n "; LATEX.write(string); 
string = "A plot of the immersed boundary at time $t=0$ and $t = " + str(T) + "$ is given in Figure \\ref{Fig:FirstAndLast}."; LATEX.write(string);
string = "\\begin{figure}[h] \n"; LATEX.write(string); 
string = "\\begin{center} \n"; LATEX.write(string);  
string = "\\includegraphics[width=2.65in, height=2.65in, angle=0]{" +str(FirstAndLastName) + "} \\\\ \n"; LATEX.write(string)
string = "\\caption{Plot of immersed boundary at $t = 0$ and $t = " + str(T) + "$. \\label{Fig:FirstAndLast}} \n"; LATEX.write(string)
string = "\\end{center} \n"; LATEX.write(string)
string = "\\end{figure} \n"; LATEX.write(string)
# Plot of the FFTs of different Modes. 
string = "\\section{Plot of FFT Modes: }\n\n "; LATEX.write(string); 
string = "Equally distributing a set of points about the fiber with respect to $\\theta$, the amplitude of the " + str(int(bub_p))
string = string + " mode is tracked and plotted vs. time and shown as a red line in Figure \\ref{Fig:FFTModep}."; LATEX.write(string);
string = "\\begin{figure}[h] \n"; LATEX.write(string); 
string = "\\begin{center} \n"; LATEX.write(string);  
string = "\\includegraphics[width=6in, height=2in, angle=0]{AmpModePlot3.eps} \\\\ \n"; LATEX.write(string)
string = "\\caption{Plot of the FFT of mode " + str(int(bub_p)) + ". \\label{Fig:FFTModep}} \n"; LATEX.write(string)
string = "\\end{center} \n"; LATEX.write(string)
string = "\\end{figure} \n"; LATEX.write(string)
string = "The best fit parameters $d_r$ and $\\omega$ are found to fit the tracted FFT mode such that \n\\[ y = e^{(-d_r t)} \\cos(\\omega t). \\]\n "; LATEX.write(string);
string = "Specifically the code found \n \\[ d_r = " + str(rate) + " \\hskip0.5in \\mbox{and } \\hskip0.5in \\omega = " + str(frq) + ".\\] \n "; LATEX.write(string); 
string = "The best fit function is plotted as a green line in Figure \\ref{Fig:FFTModep}. "; LATEX.write(string);   
string = "A plot of a second excited mode is also tracked and plotted. Specifically the " + str(int(2*bub_p)) + " mode is shown in Figure \\ref{Fig:FFT2Modep}. "; LATEX.write(string); 
string = "\\begin{figure}[h] \n"; LATEX.write(string); 
string = "\\begin{center} \n"; LATEX.write(string);  
string = "\\includegraphics[width=6in, height=2in, angle=0]{Amp2ModePlot3.eps} \\\\ \n"; LATEX.write(string)
string = "\\caption{Plot of the FFT of mode " + str(2*int(bub_p)) + ". \\label{Fig:FFT2Modep}} \n"; LATEX.write(string)
string = "\\end{center} \n"; LATEX.write(string)
string = "\\end{figure} \n"; LATEX.write(string)

# Section Report Location
string = "\\section{Report Location:} \n"; LATEX.write(string); 
string = "\n\n" + "This report was generated using the data located in: \\\\ \n{\\footnotesize \n \\verb' " + str(DirectoryShort) + " '}\n"; LATEX.write(string)
LATEX.write("\\end{document} \n")
LATEX.close()


# =========================================================================
# Write FFTData.m file for further use of the found FFT coefficents 
# =========================================================================
matlabfile = DirectoryOut + matlabfile
MATLAB = open(matlabfile,"w")
MATLAB.write("% FFT values for the given mode and 2 times the given mode as well as times. \n") 
MATLAB.write("\n\n")
# write out times
string = "dt = " + str(dt) + "; \n"
MATLAB.write(string)
MATLAB.write("t = [ ")
for dataset in range(size(amp_pmode[:,0])):
    string = str(time[dataset,0]) + ';' 
    MATLAB.write(string)
MATLAB.write("]; \n\n")
# write out FFT pmode
string1 = "pmode_Re" + str(Re) + "lambda" + str(We) + "pmode" + str(bub_p) + " = [";
MATLAB.write(string1);
for dataset in range(size(amp_pmode[:,0])):
    string = str(amp_pmode[dataset,0]) + ';'  
    MATLAB.write(string)
MATLAB.write("]; \n\n")
# write out FFT 2pmode
string2 = "pmode2_Re" + str(Re) + "lambda" + str(We) + "pmode" + str(bub_p) + " = [";
MATLAB.write(string2);
for dataset in range(size(amp_pmode[:,0])):
    string = str(amp_2pmode[dataset,0]) + ';'  
    MATLAB.write(string)
MATLAB.write("]; \n\n")
string = "plot(t," + str(string1) + ");\n"; 
MATLAB.write(string)
MATLAB.write("figure\n")
string = "plot(t," + str(string2) + ");\n";
MATLAB.write(string)


