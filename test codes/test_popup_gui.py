from tkinter import scrolledtext
import sys
from tkinter import filedialog
import tkinter as tk


########################################################################
class RedirectText(object):
    """"""

    # ----------------------------------------------------------------------
    def __init__(self, text_ctrl):
        """Constructor"""
        self.output = text_ctrl

    # ----------------------------------------------------------------------
    def write(self, string):
        """"""
        new_window = tk.Toplevel()
        label = tk.Label(new_window, text=string)
        label.pack()
        #  self.output.insert(tk.END, string)

    def flush(self):
        pass


########################################################################
class MyApp(object):
    """"""

    # ----------------------------------------------------------------------
    def __init__(self, parent):
        """Constructor"""
        self.root = parent
        self.root.title("Redirect")
        self.frame = tk.Frame(parent)
        self.frame.pack()

        self.text = scrolledtext.ScrolledText(self.frame)
        self.text.pack()

        # redirect stdout
        redir = RedirectText(self.text)
        sys.stdout = redir

        btn = tk.Button(self.frame, text="Open file", command=self.open_file)
        btn.pack()

        btn2 = tk.Button(self.frame, text="Open file", command=self.t)
        btn2.pack()

    # ----------------------------------------------------------------------
    def open_file(self):
        """
        Open a file, read it line-by-line and print out each line to
        the text control widget
        """
        options = {}
        options['defaultextension'] = '.txt'
        options['filetypes'] = [('all files', '.*'), ('text files', '.txt')]
        options['initialdir'] = '/home'
        options['parent'] = self.root
        options['title'] = "Open a file"

        with filedialog.askopenfile(mode='r', **options) as f_handle:
            for line in f_handle:
                print(line)

    def t(self):
        diciontary = {'key': value}
        diciontary.pop('test')

# ----------------------------------------------------------------------
if __name__ == "__main__":
    root = tk.Tk()
    root.geometry("800x600")
    app = MyApp(root)
    root.mainloop()