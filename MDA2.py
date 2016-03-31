import matplotlib.pyplot as plt
from math import ceil
from astroML.stats import sigmaG
import numpy as np
from uncertainties import ufloat
from astroML.stats import binned_statistic_2d
from collections import OrderedDict
from copy import deepcopy
from Tkinter import *
import os
import shutil
import ttk
import tkFileDialog
from astroML.plotting.mcmc import convert_to_stdev
############################################################################

master = Tk()
master.wm_title("Multi Dimensional Analyzer")
sizes=10
labele0=OrderedDict([('$\\chi^2$', 1), ('$T_{per}$', 3), ('$e$', 4), ('$\\omega$', 5), ('$K_1$', 6), ('$K_2$', 7)])
###################################
def clickAbout(toplevel):
    toplevel.wm_title("About")
    ABOUT_TEXT="Multi Dimensional Analyzer 2.0. (7/3/2016)\n Made using Numpy, PyPlot and AstroML. \n Please report bugs to kresimir.tisanic@gmail.com"
    about = Label(toplevel, text=ABOUT_TEXT, height=0, width=50)
    about.pack()
def Help():
    toplevel=Toplevel(master)
    toplevel.wm_title("How To Use?")
    ABOUT_TEXT="MDA 2.0 Takes an ASCII file as input.\n Columns should be separated by tabs or spaces. \n Upon hitting Run, it computes the statistic \n on parameters defined as columns."
    about = Label(toplevel, text=ABOUT_TEXT, height=0, width=50)
    about.pack()
def reportwin(namel,report,reportl):
    save_report=open(str(direktorij+'/report.tex'),'w')
    save_report.write(report)
    save_report.close()
    plt.figure(figsize=(4,3))
    ax=plt.gca()
    plt.axis('off')
    plt.table(cellText=reportl, colLabels=namel,loc='center')
    plt.savefig(str(direktorij+'/report.png'))
    plt.close()

def openf():
	dialog_output=tkFileDialog.askopenfilename(parent=master,title='Choose a file')
	global f
	global txtv
	global txtl
	global txtu
	global Number_of_rows
	global Number_of_columns
	global labele0
	f=open(dialog_output)
	Number_of_rows=sum([1 for sumacija in f])
	f.seek(0)
	Number_of_columns=max([len(sumacija.strip().split()) for sumacija in f])
	f.seek(0)
	separator = Frame(master,bd=1, relief=SUNKEN)
	separator.grid(row=1,column=0)
	txtv=[StringVar() for j in range(Number_of_columns)]
	txtl=[StringVar() for j in range(Number_of_columns)]
	txtu=[StringVar() for j in range(Number_of_columns)]
	Label(separator ,text='n').grid(row=1,column=0)
	Label(separator ,text='Names').grid(row=1,column=1)
	Label(separator ,text='Lower cuts').grid(row=1,column=2)
	Label(separator ,text='Upper cuts').grid(row=1,column=3)
	for k,j in labele0.items():
			if(j<Number_of_columns):			
				txtv[j].set(k)
	for i in range(Number_of_columns):
		Label(separator ,text=str(i)).grid(row=2+i,column=0)
		trt=0
		b=Entry(separator ,width=sizes,textvariable=txtv[i])
		b.grid(row=2+i,column=1)
		b=Entry(separator ,width=sizes,textvariable=txtl[i])
		b.grid(row=2+i,column=2)
		b=Entry(separator,width=sizes,textvariable=txtu[i])
		b.grid(row=2+i,column=3)

def openfolder(varijabla):
	global direktorij	
	varijabla.set(tkFileDialog.askdirectory(parent=master,title='Choose output folder'))
	direktorij=varijabla.get()	
	if not os.path.exists(direktorij):
	   		 os.makedirs(direktorij)
def savelimits(labele,labeleu,labelel):
	savediag=Toplevel(master)
	save_name=StringVar()
	Label(savediag,text="Save file:").grid(row=0,column=0)
	Entry(savediag,width=sizes,textvariable=save_name).grid(row=0,column=1)
	Button(savediag,text='Save',command=lambda:start_save(save_name.get())).grid(row=0,column=2)
	def start_save(save_name):
		filesave=open(str(save_name+'.dat'),'w')
		for i in range(len(labele)):
			filesave.write(str(str(i)+':'+str(labele[i])+':'+str(labelel[i])+':'+str(labeleu[i])+'\n'))
		filesave.close()
		savediag.destroy()
def open_limits():	
	opendiag=Toplevel(master)
	save_name=StringVar()
	Label(opendiag,text="Open save file:").grid(row=0,column=0)
	Entry(opendiag,width=sizes,textvariable=save_name).grid(row=0,column=1)
	Button(opendiag,text='Open',command=lambda:start_open(save_name)).grid(row=0,column=2)
	def start_open(save_name):
		filesave=open(str(str(save_name.get())+'.dat'),'r')
		global txtv
		global txtu
		global txtl
		j=0
		for line in filesave:
			line1=line.strip()
			line2=line1.split(':')
			txtv[j].set(line2[1])
			txtl[j].set(line2[2])
			txtu[j].set(line2[3])
			j=j+1
		filesave.close()
		opendiag.destroy()

def sigmas(axis1, mean1,sigma1,axis2,mean2,sigma2):
		    return np.sqrt(((axis1-mean1)/sigma1)**2.+((axis2-mean2)/sigma2)**2.)
def start():
	direktorij=output_folder.get()	
	if not os.path.exists(direktorij):
	   		 os.makedirs(direktorij)
	finishing_string='Name & Mean & Median\\\\ \n'
	starting_list=['Name','Mean','Median']
	finishing_list=[]
	f.seek(0)
	textval=[a.get() for a in txtv]
	textu=[a.get() for a in txtu]
	textl=[a.get() for a in txtl]
	labele=OrderedDict([(textval[j],j) for j in range(len(textval)) if(len(textval[j])!=0)])
	labeleu=OrderedDict([(textval[j],float(textu[j])) for j in range(len(textval)) if((len(textu[j])!=0) and (len(textval[j])!=0))])
	labelel=OrderedDict([(textval[j],float(textl[j])) for j in range(len(textval)) if((len(textl[j])!=0) and (len(textval[j])!=0))])	
	X=np.array([[0. for i in range(Number_of_columns)] for j in range(Number_of_rows)])
	i=0
	for line in f:
			line0=line.strip()
			line1=line0.split()
			passing=1
			for dicts, vals in labele.items():
				if dicts in labelel:
					if(labelel[dicts]>float(line1[vals])):
						passing=0
				if dicts in labeleu:
					if(labeleu[dicts]<float(line1[vals])):
						passing=0
			if(passing==1):			
				for j in range(Number_of_columns):
					X[i][j]=float(line1[j])
				
				i=i+1
	X=X[0:i,:]
	X=np.array([[X[k][l] for l in range(Number_of_columns)] for k in range(i) ])
	fig=plt.figure(figsize=(20,15))
	count=1
	for ime, val in labele.items():
			ax=fig.add_subplot(2,3,count)
			ax.set_xlabel(ime, size=30)
			ax.hist(X[:,val], bins=50, normed=False, histtype='stepfilled',color='blue',facecolor='blue')
			ax.axvline(np.mean(X[:,val]), color='orange', linestyle='--')
			ax.axvline(np.median(X[:,val]), color='green', linestyle='--')
			ax.xaxis.major.formatter._useMathText = True
			ax.ticklabel_format(style='sci', axis='x', scilimits=(-5,5))
			text='Mean and median:\n'+str("mean$\\rightarrow$ $ {:.2uL}$\n".format(ufloat(np.mean(X[:,val]),np.std(X[:,val])))+"median$\\rightarrow$ $ {:.2uL}$".format(ufloat(np.median(X[:,val]),sigmaG(X[:,val]) ) ))
			finishing_string=finishing_string+ime+" & "+"$ {:.2uL}$".format(ufloat(np.mean(X[:,val]),np.std(X[:,val])))+" & $ {:.2uL}$".format(ufloat(np.median(X[:,val]),sigmaG(X[:,val])))+"\\\\ \n"			
			finishing_list.append([ime,"$ {:.2uL}$".format(ufloat(np.mean(X[:,val]),np.std(X[:,val])))," $ {:.2uL}$".format(ufloat(np.median(X[:,val]),sigmaG(X[:,val])))])
			ax.text(.65,.9,text,transform = ax.transAxes)
			count=count+1
	plt.tight_layout()
	plt.savefig(direktorij+'/Histogrami.png')
	plt.close()
	reportwin(starting_list,finishing_string,finishing_list)
	for names0,vals0 in labele.items():
				fig=plt.figure(figsize=(30,25))
				labele1=deepcopy(labele)

				del labele1[names0]
				labele2=deepcopy(labele1)
				nx=ceil(np.sqrt(len(labele1)*(len(labele1)-1)/2))
				ny=ceil(len(labele1)*(len(labele1)-1)/2/nx)
				#fig, axes = plt.subplots(nrows=2, ncols=2, sharex=True, sharey=True)
				counts=1
				cmap_multicolor = plt.cm.jet
				
				for names,vals in labele1.items():
					del labele2[names]
					for names1, vals1 in labele2.items():
						N0, xedges0, yedges0 = binned_statistic_2d(X[:,vals], X[:,vals1], X[:,labele[names0]], 'mean', bins=100)
						ax=fig.add_subplot(ny,nx,counts)
						im=ax.imshow(N0.T, origin='lower',extent=[xedges0[0], xedges0[-1], yedges0[0], yedges0[-1]], aspect='auto', interpolation='nearest', cmap=cmap_multicolor)
						plt.xlim(xedges0[0], xedges0[-1])
						plt.ylim(yedges0[0], yedges0[-1])
						plt.xlabel(names, size=30)
						plt.ylabel(names1, size=30)
						ax.xaxis.major.formatter._useMathText = True
						ax.yaxis.major.formatter._useMathText = True
						ax.ticklabel_format(style='sci', axis='x', scilimits=(-5,5))
						#m_1 = np.linspace(xedges0[0], xedges0[-1], 100)
						#m_2 = np.linspace(yedges0[0], yedges0[-1], 100)

						#MX,MY = np.meshgrid(m_1, m_2)

						#Z = sigmas(MX,np.median(X[:,vals]),sigmaG(X[:,vals]), MY,np.median(X[:,vals1]),sigmaG(X[:,vals1]))

						H, xbins, ybins = np.histogram2d(X[:,vals], X[:,vals1],bins=100)

						Nsigma = convert_to_stdev(np.log(H))
						cont=plt.contour(0.5 * (xbins[1:] + xbins[:-1]),0.5 * (ybins[1:] + ybins[:-1]),Nsigma.T,levels=[0.6827,0.6827,0.9545, 0.9545], colors=['.25','.25','0.5','0.5'],linewidths=2)					
						counts=counts+1
				
				cmap_multicolor.set_bad('w', 1.)
				fig.subplots_adjust(bottom=0.1)
				cbar_ax = fig.add_axes([0.1, 0.05, 0.8, 0.025])
				cb=fig.colorbar(im, cax=cbar_ax, format=r'$%.1f$',orientation='horizontal')
				cb.set_label(str('$\\langle '+names0.replace('$','')+'\\rangle $'), size=30)
				
				plt.savefig(direktorij+'/'+''.join([i for i in names0 if (i.isalpha() or i.isdigit())])+'.png',bbox_inches='tight')
				plt.close()

output_folder=StringVar(value='results')


   
menubar = Menu(master)
filemenu = Menu(menubar, tearoff=0)
filemenu.add_command(label="Open",command=openf)
filemenu.add_command(label="Output folder",command=lambda:openfolder(output_folder))
filemenu.add_separator()
filemenu.add_command(label='Open limits', command=lambda:open_limits())
filemenu.add_command(label='Save limits', command=lambda:savelimits([a.get() for a in txtv],[a.get() for a in txtu],[a.get() for a in txtl]))
filemenu.add_separator()
filemenu.add_command(label='Exit', command=lambda:master.quit())

helpmenu=Menu(menubar, tearoff=0)
helpmenu.add_command(label="About",command=lambda:clickAbout(Tk()))
helpmenu.add_command(label="How to use?",command=Help)
menubar.add_cascade(label="File", menu=filemenu)
menubar.add_command(label="Run",command=lambda:start())
menubar.add_cascade(label="Help", menu=helpmenu)
master.config(menu=menubar)


mainloop()
