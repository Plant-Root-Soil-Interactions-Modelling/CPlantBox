from tkinter import *
from scipy.special.tests.test_data import eval_genlaguerre_ddd


def NewFile():
    print("New File!")


def OpenFile(): 
    print("na")


def About():
    print("This is a simple example of a menu")


def Quit():
    print("quite")

    
root = Tk()

s = tkk.Style(root)

menu = Menu(root)
root.config(menu=menu)
filemenu = Menu(menu)
menu.add_cascade(label="File", menu=filemenu)
filemenu.add_command(label="New", command=NewFile)
filemenu.add_command(label="Open...", command=OpenFile)
filemenu.add_separator()
filemenu.add_command(label="Exit", command=root.quit)

helpmenu = Menu(menu)
menu.add_cascade(label="Help", menu=helpmenu)
helpmenu.add_command(label="About...", command=About)

button = Button(master=root, text="Quit", command=Quit())
button.pack(side=BOTTOM)

mainloop()
