#!/usr/bin/env python
from tkinter import ttk
import os
from tkinter import *
import re
import time
from tkinter import filedialog as fd
import mysql.connector 

class InputForm():
	def __init__ (self):
		self.filelist3=[]
		self.filelist=""
		self.PO=""
		self.kit=""		

		self.main = Tk()
		self.main.title("RNAseq_Pipline")
		self.main.geometry("1000x800")
		
		self.frame0 = ttk.Frame(self.main, padding=(3, 3, 12, 12))#
		self.frame0.grid(column=0, row=0, padx=10, pady=2,columnspan=2, sticky=(N, S, E, W))#, sticky=(N, S, E, W)
		
		Label(self.frame0,text = "Path:").grid(row = 0, column = 0, padx=70, sticky=W)
		self.path = StringVar()
		Entry(self.frame0, textvariable = self.path, width=50).grid(row = 0, column = 1, sticky=E)
		Button(self.frame0, text = "select", command = self.callback).grid(row = 0, column = 2, padx=10, sticky=W)

		Label(self.frame0,text = "Species").grid(row = 1, column = 0, pady=4, padx=12, sticky=(N, S, E, W))#
		Label(self.frame0,text = "library").grid(row = 1, column = 1, pady=4, padx=12, sticky=(N, S, E, W))#
		Label(self.frame0,text = "base").grid(row = 1, column = 2, pady=4, padx=12, sticky=(N, S, E, W))#
		
		self.path_0 = StringVar()
		self.path_1 = StringVar()
		self.path_2 = StringVar()
		Entry(self.frame0, textvariable = self.path_0, width=10).grid(row = 2, column = 0, pady=4, padx=12, sticky=(N, S, E, W))#
		Entry(self.frame0, textvariable = self.path_1, width=10).grid(row = 2, column = 1, pady=4, padx=12, sticky=(N, S, E, W))#
		Entry(self.frame0, textvariable = self.path_2, width=10).grid(row = 2, column = 2, pady=4, padx=12, sticky=(N, S, E, W))#
		Button(self.frame0, text = "information", command = self.callback2).grid(row = 2, column = 3, pady=4, padx=12, sticky=(N, S, E, W))#
		
		self.frame = ttk.Frame(self.main, padding=(3, 3, 12, 12))#
		self.frame.grid(column=0, row=1, padx=10, pady=2, sticky=(N, S, E, W))#, sticky=(N, S, E, W)
		self.frame2 = ttk.Frame(self.main, padding=(3, 3, 12, 12))
		self.frame2.grid(column=1, row=1, padx=10, pady=2, sticky=(N, S, E, W))#, padx=10, pady=2, sticky=(N)
		
		self.label1 =ttk.Treeview(self.frame2,columns=['1','2'],show='headings')
		self.label1.heading('1',text='Sample')
		self.label1.heading('2',text='Group')
		self.label1.grid()
		
		Label(self.frame, text='Sample list', width=50).grid(column=0, row=0, columnspan=2)
		valores = StringVar()
		self.lstbox = Listbox(self.frame, listvariable=valores, selectmode=MULTIPLE, width=50, height=10)
		self.lstbox.grid(column=0, row=1, columnspan=2)					
		Label(self.frame, text='Group name').grid(column=0, row=2, columnspan=2)
		self.var=StringVar()
		Entry(self.frame,textvariable=self.var, width=50).grid(column=0, row=3, columnspan=2)
		self.btn = Button(self.frame, text="Group", command=self.select, width=25)
		self.btn.grid(column=0, row=4, pady=4)	
		self.btn5 = Button(self.frame, text="clear_group", command=self.cleardata, width=25)
		self.btn5.grid(column=1, row=4, pady=4, padx=12)	
		
		self.frame3 = ttk.Frame(self.main, padding=(3, 3, 12, 12))
		self.frame3.grid(column=0, row=2, padx=10, pady=2, sticky=(N, S, E, W))
		self.frame5 = ttk.Frame(self.main, padding=(3, 3, 12, 12))
		self.frame5.grid(column=1, row=2, padx=10, pady=2, sticky=(N, S, E, W))
		
		self.label2 =ttk.Treeview(self.frame5,columns=['1','2'],show='headings')
		self.label2.heading('1',text='Treatment')
		self.label2.heading('2',text='control')
		self.label2.grid()
		
		Label(self.frame3, text='Treatment').grid(column=0, row=0)
		self.lstbox2 = ttk.Combobox(self.frame3, values=self.filelist3, width=23, height=10)				
		self.lstbox2.grid(column=0, row=1, pady=4)
		
		Label(self.frame3, text='Control').grid(column=1, row=0)
		self.lstbox3 = ttk.Combobox(self.frame3, values=self.filelist3, width=23, height=10)				
		self.lstbox3.grid(column=1, row=1, pady=4, padx=12)
		
		self.btn2 = Button(self.frame3, text="compare", command=self.select2, width=25)
		self.btn2.grid(column=0, row=2, pady=4)
		
		self.btn3 = Button(self.frame3, text="clear_compare", command=self.cleardata_2, width=25)
		self.btn3.grid(column=1, row=2, pady=4, padx=12)
		
		self.btn4 = Button(self.frame3, text="output", command=self.savedata, width=54)
		self.btn4.grid(column=0, row=3, padx=10, pady=2,columnspan=2)#
		
		Label(self.frame3, text='Species').grid(column=0, row=4)
		self.combo0 = ttk.Combobox(self.frame3, values=["1.H.sapiens","2.M.musculus","3.R.norvegicus","4.A.thaliana","5.C.elgean","6.D.rerio","7.D.melanogaster","8.O.sativa","9.Custom"],state="readonly", width=23, height=10)
		self.combo0.grid(column=0, row=5)
		self.combo0.current(0)
		
		Label(self.frame3, text='FC').grid(column=1, row=4, sticky=E)
		self.FC=DoubleVar()
		self.FC.set(2)
		Entry(self.frame3,textvariable=self.FC, width=13).grid(column=1, row=5, sticky=E)
		
		Label(self.frame3, text='P-value').grid(column=1, row=4, sticky=W)
		self.PV=DoubleVar()
		self.PV.set(0.05)
		Entry(self.frame3,textvariable=self.PV, width=13).grid(column=1, row=5, sticky=W)
		
		Label(self.frame3, text='Lib-strandness').grid(column=0, row=6)
		self.combo = ttk.Combobox(self.frame3, values=["1.fr-firststrand","2.fr-secondstrand","3.fr-unstranded"],state="readonly", width=23, height=10)
		self.combo.grid(column=0, row=7)
		self.combo.current(0)

		Label(self.frame3, text='normalization').grid(column=1, row=6)
		self.combo1 = ttk.Combobox(self.frame3, values=["0.TPM","1.FPKM"],state="readonly", width=23, height=10)
		self.combo1.grid(column=1, row=7)
		self.combo1.current(0)
		
		Label(self.frame3, text='Other analysis').grid(column=0, row=8, columnspan=2)
		self.Fusion = BooleanVar()
		Checkbutton(self.frame3,text = "STAR-Fusion",variable = self.Fusion).grid(column=0, row=9, sticky = W)
		self.AS = BooleanVar()
		Checkbutton(self.frame3,text = "Alternative Splicinng(MISO)",variable = self.AS).grid(column=1, row=9, sticky = W)

		self.Fusion_FFPE = BooleanVar()
		Checkbutton(self.frame3,text = "STAR-Fusion_FFPE",variable = self.Fusion_FFPE).grid(column=0, row=10, sticky = W)
		self.arriba_FFPE = BooleanVar()
		Checkbutton(self.frame3,text = "arriba_FFPE",variable = self.arriba_FFPE).grid(column=1, row=10, sticky = W)

		self.MATS = BooleanVar()
		Checkbutton(self.frame3,text = "MATS",variable = self.MATS).grid(column=0, row=11, sticky = W)
		self.other = BooleanVar()
		Checkbutton(self.frame3,text = "other",variable = self.other).grid(column=1, row=11, sticky = W)
		
		self.btn6 = Button(self.frame3, text="output-RNAseq common line", command=self.savedata2, width=54)
		self.btn6.grid(column=0, row=12, padx=12, pady=2,columnspan=2)
	
		self.frame6 = ttk.Frame(self.main, padding=(3, 3, 12, 12))
		self.frame6.grid(column=0, row=3, padx=10, pady=2,columnspan=2, sticky=(N, S, E, W))
		self.Output = Text(self.frame6,width=130, height=8,bg = "light cyan")
		self.Output.grid(column=0, row=0,columnspan=2)
	
		self.main.mainloop()
	
	def callback2(self):

		m=self.PO
		db=mysql.connector.connect(database='LIMS',host='192.168.5.114',user='ngs',passwd='bio@NGS80158777')
		mycursor=db.cursor()
		
		print(m)
		
		mycursor = db.cursor()
		mycursor.execute("Select * from `POrecord` where `poNUM`=\'"+m+"\'")
		myresult0=mycursor.fetchall()
	
		mycursor = db.cursor()
		mycursor.execute("Select * from `projects` where `id`=\'"+str(myresult0[0][4])+"\'")
		myresult2=mycursor.fetchall()


		mycursor = db.cursor()
		mycursor.execute("Select * from `library` where `POrecord_id`=\'"+str(myresult0[0][0])+"\'")
		myresult=mycursor.fetchall()

		mycursor = db.cursor()
		mycursor.execute("Select * from `sequencing` where `POrecord_id`=\'"+str(myresult0[0][0])+"\'")
		myresult1=mycursor.fetchall()

		self.path_0.set(myresult2[0][7])
		self.path_1.set(myresult[0][3])
		self.path_2.set(str(myresult1[0][10])+str(myresult1[0][11]))
		self.kit=myresult[0][3]
	
	def callback(self):
		self.lstbox.delete(0,END)
		text =fd.askdirectory()
		self.filelist=text
		self.filelist3=[]
		text1=os.path.dirname(text)
		self.PO=os.path.basename(text1)
		#reg=re.compile(r'NGS[0-9]+')
		#matc=reg.search(text1)
		print(self.PO)

		for root, dirs, files in os.walk(text1, topdown=False):
			for name in files:
				m=re.match(r"(.+)_S[0-9]+_?L?[0-9]*_R1_[0-9]+.fastq.gz",name)
				
				str=os.path.join(root, name)
				if m:		
					self.filelist3.append(m.group(1))
				if str.split('.')[-1]=='bam':
					str1=os.path.split(str)[1].split('.')[0]
					self.filelist3.append(str1)
		self.filelist3=set(self.filelist3)
		self.filelist3=list(self.filelist3)
		for i in self.filelist3:
			self.lstbox.insert(END, i)
		self.lstbox2["values"]=self.filelist3
		self.lstbox3["values"]=self.filelist3
		self.path.set(text)
		self.label1.delete(*self.label1.get_children())
		self.label2.delete(*self.label2.get_children())
		self.path_0.set("")
		self.path_1.set("")
		self.path_2.set("")
		self.kit=""
	def savedata2(self):
		self.Output.delete("1.0","end")
		sap=self.combo0.get().split('.')[0]
		stra=self.combo.get().split('.')[0]
		nor=self.combo1.get().split('.')[0]
		FC=str(self.FC.get())
		PV=str(self.PV.get())
		dir=os.path.dirname(self.filelist)
		text1='cd '+dir+'\n'
		if re.match(r'Others',self.kit):
			text1=text1+'rand_trim.pl 1\n'+'cd '+self.filelist+'\nln -s ../trimmedFASTQ/*.gz ./\n'
		else :
			text1=text1+'rand_trim.pl 0\n'+'cd '+self.filelist+'\nln -s ../trimmedFASTQ/*.gz ./\n'
		#text1=text1+'rand_trim.pl\n'+'cd '+self.filelist+'\nln -s ../trimmedFASTQ/*.gz ./\n'
		text1=text1+'featurecount2_G_v5.pl -s '+sap+' -l '+stra+' -f '+FC+' -p '+PV+' -e '+nor+' -r 0'
		if self.Fusion.get() | self.AS.get() |  self.MATS.get() | self.Fusion_FFPE.get() | self.arriba_FFPE.get():
			Fusion="N"
			AS="N"
			MATS="N"
			Fusion_FFPE="N"
			arriba_FFPE="N"
			if self.Fusion.get():
				Fusion="Y"
			if self.AS.get():
				AS="Y"
			if self.Fusion_FFPE.get():
				Fusion_FFPE="Y"
			if self.arriba_FFPE.get():
				arriba_FFPE="Y"
			if self.MATS.get():
				MATS="Y"
			text1=text1+' -gf '+Fusion+' -as '+AS+' -gfF '+Fusion_FFPE+' -af '+arriba_FFPE+' -rMATS '+MATS+' -sP Y'
		tt=str(time.strftime("%Y_%m_%d_%H_%M_%S", time.localtime()))
		text1=text1+' > featurecount2_G_v5_'+tt+'.log'
		self.Output.insert(END,text1)
		data=open(self.filelist+"/RNA_run"+tt+".sh",'w+')
		print(text1,file=data)
		data.close()
		
	def savedata(self):
		items = self.label1.get_children()
		if items:
			data=open(self.filelist+"/sample2GroupInfo.txt",'w+')
			for item1 in items:
				item_text = self.label1.item(item1,"values")
				print(item_text[0]+'\t'+item_text[1],file=data)
			data.close()
		data=open(self.filelist+"/groupCompareInfo.txt",'w+')
		items = self.label2.get_children()
		for item1 in items:
			item_text = self.label2.item(item1,"values")
			print(item_text[0]+'\t'+item_text[1],file=data)
		data.close()
		
	def cleardata(self):
		items = self.label1.selection()	
		for item in items:
			self.label1.delete(item)

	
	def cleardata_2(self):
		items = self.label2.selection()	
		for item in items:
			self.label2.delete(item)
	
	def select(self):
		reslist = list()
		seleccion = self.lstbox.curselection()
		filelist2=""
		for i in seleccion:
			entrada = self.lstbox.get(i)
			reslist.append(entrada)
		for val in reslist:
			filelist2=filelist2+val+'\t'+self.var.get()+'\n'
			li=[val,self.var.get()]
			self.label1.insert('','end',values=li)
		self.filelist3.append(self.var.get())
		self.filelist3=set(self.filelist3)
		self.filelist3=list(self.filelist3)
		self.lstbox2["values"]=self.filelist3
		self.lstbox3["values"]=self.filelist3
	
	def select2(self):
		reslist = list()
		seleccion1 = self.lstbox2.get()
		seleccion2 = self.lstbox3.get()
		li=[seleccion1,seleccion2]
		self.label2.insert('','end',values=li)

abc =InputForm().filelist

print("returned value is:", abc)
