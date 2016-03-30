from Tkinter import *
from math import *
import sys
import urllib
import urllib2
from numpy import sign
import tkMessageBox
import smtplib
import webbrowser
from threading import Timer

class App:
    sizes=10
    hoursconvert=15.
    radiansconvert=pi/180
    def __init__(self,master):
        # Observational parameters
        self.horizontal=LabelFrame(master,text="Observation parameters",width=300,height=120)
        self.horizontal.grid(row=0,column=0)
        self.horizontal.grid_propagate(False)
        
        Label(self.horizontal,text='Lattitude').grid(row=1,column=0)
        Label(self.horizontal,text='deg/h').grid(row=0,column=1)
        Label(self.horizontal,text='min').grid(row=0,column=2)
        Label(self.horizontal,text='sec').grid(row=0,column=3)
        
        self.h1 = Entry(self.horizontal,width=self.sizes)
        self.h1.grid(row=1,column=1)
        self.h1.insert(0, '0')
        self.h2=Entry(self.horizontal,width=self.sizes)
        self.h2.grid(row=1,column=2)
        self.h2.insert(0, '0')
        self.h3=Entry(self.horizontal,width=self.sizes)
        self.h3.grid(row=1,column=3)
        self.h3.insert(0, '0')
        Label(self.horizontal,text='Longitude').grid(row=2,column=0)
        
        self.h4 = Entry(self.horizontal,width=self.sizes)
        self.h4.grid(row=2,column=1)
        self.h4.insert(0, '0')
        self.h5=Entry(self.horizontal,width=self.sizes)
        self.h5.grid(row=2,column=2)
        self.h5.insert(0, '0')
        self.h6=Entry(self.horizontal,width=self.sizes)
        self.h6.grid(row=2,column=3)
        self.h6.insert(0, '0')
        
        Button(self.horizontal,text='Sidereal Time',command=self.fetch).grid(row=3,column=0)
        self.st1=Entry(self.horizontal,width=self.sizes)
        self.st2=Entry(self.horizontal,width=self.sizes)
        self.st3=Entry(self.horizontal,width=self.sizes)
        self.st1.insert(0,'0')
        self.st2.insert(0,'0')
        self.st3.insert(0,'0')
        self.st1.grid(row=3,column=1)
        self.st2.grid(row=3,column=2)
        self.st3.grid(row=3,column=3)
        
        # Input RA/DEC
        self.frame=LabelFrame(master,text="Right ascension and declination",width=300,height=90)
        self.frame.grid(row=1,column=0)
        self.frame.grid_propagate(False)
        
        Label(self.frame,text='deg/h').grid(row=0,column=1)
        Label(self.frame,text='min').grid(row=0,column=2)
        Label(self.frame,text='sec').grid(row=0,column=3)
        
        
        Label(self.frame,text='RA').grid(row=1,column=0)
        self.e = Entry(self.frame,width=self.sizes)
        self.e.grid(row=1,column=1)
        self.e.insert(0, '0')
        self.e1=Entry(self.frame,width=self.sizes)
        self.e1.grid(row=1,column=2)
        self.e1.insert(0, '0')
        self.e2=Entry(self.frame,width=self.sizes)
        self.e2.grid(row=1,column=3)
        self.e2.insert(0, '0')
        
        Label(self.frame,text='DEC').grid(row=2,column=0)
        self.e3 = Entry(self.frame,width=self.sizes)
        self.e3.grid(row=2,column=1)
        self.e3.insert(0, '0')
        self.e4=Entry(self.frame,width=self.sizes)
        self.e4.grid(row=2,column=2)
        self.e4.insert(0, '0')
        self.e5=Entry(self.frame,width=self.sizes)
        self.e5.grid(row=2,column=3)
        self.e5.insert(0, '0')
        
        ### Output
        self.output=LabelFrame(master,width=300,height=160)
        self.output.grid(row=2,column=0)
        self.output.grid_propagate(False)
        
        #Calculate, radians and bug report buttons
        self.var = IntVar()
        lower_frame=LabelFrame(master,width=300,height=30)
        lower_frame.grid(row=3,column=0)
        lower_frame.grid_propagate(False)
        Checkbutton(lower_frame, text="radians", variable=self.var).grid(row=0,column=0)
        Button(lower_frame, text="calculate",command=self.convert).grid(row=0,column=1)
        Button(lower_frame, text="Bugs",command=self.reply).grid(row=0,column=2)
        Button(lower_frame,text="Help",command=self.helper).grid(row=0,column=3)

    def helper(self):
        helping=Toplevel(master)
        helpf1=LabelFrame(helping,text="Input")
        helpf1.pack()
        helpf2=LabelFrame(helping,text="Output")
        helpf2.pack()
        helpf3=LabelFrame(helping,text="About")
        helpf3.pack()
        
        helping.wm_title("Help")
        helps1 = Label(helpf1, text="Enter lattitude and longitude of your observing site\n Enter siderial time in h/m/s or calculate it from website\n http://tycho.usno.navy.mil/cgi-bin/sidereal-post.sh doesn't always respond\n Enter right ascension as h/m/s\n Enter declination as deg/m/s\n Choose wheter to output in radians instead of h/deg/m/s.", height=0, width=70)
        helps1.pack()
        helps2 = Label(helpf2, text="RA=Right ascension\n DEC=declination\n z=Zenith distance\n h=altitude\n H=hour angle\n A=azimuth", height=0, width=70)
        helps2.pack()
        helps3 = Label(helpf3, text="March, 2016\n If you find bugs, please report them to kresimir.tisanic@gmail.com\n", height=0, width=70)
        helps3.pack()
        
    def reply(self):
            #send bug report
            webbrowser.open('mailto:kresimir.tisanic@gmail.com?Subject=Bug report - Astro converter', new=1)
            
    def fetch(self):
            #getting siderial time from http://tycho.usno.navy.mil/cgi-bin/sidereal-post.sh
            url = 'http://tycho.usno.navy.mil/cgi-bin/sidereal-post.sh'
            
            ### preparing data for sending            
            def signs(x):
                if sign(x)>0:
                    return "+1"
                else:
                    return "-1"
            values = {'lond' : str(abs(float(self.h4.get()))),
                      'lonm' : str(self.h5.get()),
                      'lons' : str(round(float(self.h6.get()))),
                      'ew':signs(float(self.h4.get()))}
            data = urllib.urlencode(values)

            ### Requesting response
            req = urllib2.Request(url, data)
            response = urllib2.urlopen(req)
            t = Timer(5.0, self.fetch,[response])
            t.start()
            the_page = response.read()
            t.cancel()
            y=the_page.split('LST')
            z=y[0].split("<H1>")
            
            ### outputting siderial time
            self.sidh,self.sidm, self.sids= z[1].split(":")
            self.st1.delete(0,END)
            self.st2.delete(0,END)
            self.st3.delete(0,END)            
            self.st1.insert(0,self.sidh)
            self.st2.insert(0,self.sidm)
            self.st3.insert(0,self.sids)

    def convert(self):
       self.sidh=self.st1.get()
       self.sidm=self.st2.get()
       self.sids=self.st3.get()
       #conversion function  with output
       def conversion(x,y,z):
           # converts degrees/hours to decimal form
           return (abs(float(x))+float(y)/60+float(z)/3600)*sign(float(x))
       def conversion_return(x,ch):
           # converts decimal form to degrees 
           sgn=sign(x)
           x=abs(x)
           x1=x-x%1
           x2=floor((3600*(x%1))/60)
           x3=(3600*(x%1))-60*floor((3600*(x%1))/60)
           return str( "%.0f"%(sgn*x1)+ch+" %.0f"%( x2)+'min '+"%.2f"%(x3)+'s')
       ### input conversion
       self.izrazra=conversion(self.e.get(),self.e1.get(),self.e2.get())
       self.izrazdec=conversion(self.e3.get(),self.e4.get(),self.e5.get())
       self.izrazra_radians=self.izrazra*App.radiansconvert*App.hoursconvert
       self.izrazdec_radians=self.izrazdec*App.radiansconvert
       self.phi_radians=conversion(self.h1.get(),self.h2.get(),self.h3.get())*App.radiansconvert
       self.siderial_radians=conversion(self.sidh,self.sidm,self.sids)*App.radiansconvert*App.hoursconvert
       
       ### further calculations
       self.hour_radians=self.siderial_radians-self.izrazra_radians
       self.zenithdist_radians=acos(cos(self.hour_radians)*cos(self.izrazdec_radians)*cos(self.phi_radians)+sin(self.izrazdec_radians)*sin(self.phi_radians))
       self.height_radians=pi*0.5-self.zenithdist_radians
       
       ### checking for undefined azimuth
       if(abs(sin(self.zenithdist_radians))>10**(-10) and abs(sin(self.hour_radians)*cos(self.izrazdec_radians)/sin(self.zenithdist_radians))<1.):
           self.azimuth_radians=asin(sin(self.hour_radians)*cos(self.izrazdec_radians)/sin(self.zenithdist_radians))
           if self.hour_radians<0:
               self.azimuth_radians=pi+self.azimuth_radians
               self.hour_radians=2*pi+self.hour_radians
           else:    
               self.azimuth_radians=pi-self.azimuth_radians
       else:
           self.azimuth_radians='error'
       
       ### changing behaviour if radians output is selected
       if(self.var.get()==1):
           self.output.config(text="Results in radians")
           self.izrazra=str(self.izrazra_radians)
           self.izradec=str(self.izrazdec_radians)
           self.zenithdist=str(self.zenithdist_radians)
           self.hour=str(self.hour_radians)
           self.height=str(self.height_radians)
           if(self.azimuth_radians!='error'):
               self.azimuth=str(self.azimuth_radians)
       else:
           self.output.config(text="Results in degrees")
           self.zenithdist=conversion_return(self.zenithdist_radians/App.radiansconvert,'deg')
           self.hour=conversion_return(self.hour_radians/App.radiansconvert/App.hoursconvert,'h')
           self.height=conversion_return(self.height_radians/App.radiansconvert,'deg')
           self.izrazra=conversion_return(self.izrazra,'deg')
           self.izrazdec=conversion_return(self.izrazdec,'deg')
           if(self.azimuth_radians!='error'):
               self.azimuth=conversion_return(self.azimuth_radians/App.radiansconvert,'deg')
       
       # Showing results
       w=Entry(self.output,width=4*self.sizes)
       Label(self.output,text='RA').grid(row=0,column=1)
       w.insert(0,self.izrazra)
       w.grid(row=0,column=2)
       
       w1=Entry(self.output,width=4*self.sizes)
       w1.insert(0,self.izrazdec)
       w1.grid(row=1,column=2)
       Label(self.output,text='DEC').grid(row=1,column=1)
       
       w2=Entry(self.output,width=4*self.sizes)
       w2.insert(0,self.zenithdist)
       w2.grid(row=2,column=2)
       Label(self.output,text='z').grid(row=2,column=1)
       
       w5=Entry(self.output,width=4*self.sizes)
       w5.insert(0,self.height)
       w5.grid(row=3,column=2)
       Label(self.output,text='h').grid(row=3,column=1)
       
       w4=Entry(self.output,width=4*self.sizes)
       w4.insert(0,self.hour)
       w4.grid(row=4,column=2)
       Label(self.output,text='H').grid(row=4,column=1)
       
       w3=Entry(self.output,width=4*self.sizes)
       if(self.azimuth_radians!='error'):
           w3.insert(0,self.azimuth)
       else:
           w3.insert(0,'undefined')
       w3.grid(row=5,column=2)
       Label(self.output,text='A').grid(row=5,column=1)

master = Tk()
master.wm_title("Astro converter")
app=App(master)
master.mainloop()