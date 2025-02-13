from tkinter import *
#from tkinter import messagebox
#from tkMessageBox import *
from tkinter.colorchooser import askcolor              
from tkinter.filedialog   import askopenfilename

#root = Tk()

def GetFileName(filext="*"):
    """ Returns a file name - a file has to be chosen
        Input: filext = "*" - File extention in the first list
    """
    root = Tk()
    file = ""
    if filext == "*":
        fileclass = "All Files"
        fitype = "*.*"
    if filext == "xls":
        fileclass = 'Excel File'
        fitype = "*." + filext
    while file == "":
        file = askopenfilename(filetypes=((fileclass,fitype),
                                          ("Text File","*.txt*"),
                                          ("LP Files","*.lp")),
                               title= "Choose a file")
        root.withdraw()
        print(file)      
    return file    


def ViewFileName(filext="*"):
    """ View files - Choose one or cancel
        Input: filext = "*" - File extention in the first list
    """
    root = Tk()
    file = ""
    if filext == "*":
        fileclass = "All Files"
        fitype = "*.*"
    elif filext == "xls":
        fileclass = "Excel File"
        fitype = "*." + filext
    file = askopenfilename(filetypes=((fileclass,fitype),
                                          ("Text File","*.txt*"),
                                          ("LP Files","*.lp"),
                                      ("Excel Files","*.xls")),
                               title= "Choose a file")
    root.withdraw()
    print(file)      
    return file

#file = ViewFileName()

def OpenFile():
#    from tkinter import filedialog
#    from tkinter import *

    root = Tk()
    root.filename =  filedialog.askopenfilename(title = "Select file",filetypes = (("Excel Files","*.xls"),("all files","*.*")))
    return root.filename
