import tkinter as tk
from tkinter import ttk, filedialog
from codon_analysis import aggregate_codon_counts, compute_codon_usage_from_counts, get_codon_counts_for_gene_by_locus, compare_codon_usage, validate_codon_usage_percentages
from Bio import SeqIO

class AutocompleteEntry(tk.Entry):
    def __init__(self, autocompleteList, *args, **kwargs):
        self.listbox = None
        self.autocompleteList = autocompleteList
        super().__init__(*args, **kwargs)
        self.var = self["textvariable"]
        if self.var == '':
            self.var = self["textvariable"] = tk.StringVar()

        self.var.trace('w', self.changed)
        self.bind("<Return>", self.selection)
        self.bind("<Up>", self.move_up)
        self.bind("<Down>", self.move_down)

    def changed(self, name, index, mode):
        if self.var.get() == '':
            if self.listbox:
                self.listbox.destroy()
                self.listbox = None
        else:
            words = self.comparison()
            if words:
                if not self.listbox:
                    self.listbox = tk.Listbox(width=self["width"], height=4)
                    self.listbox.bind("<Double-Button-1>", self.selection)
                    self.listbox.place(x=self.winfo_x(), y=self.winfo_y() + self.winfo_height())
                self.listbox.delete(0, tk.END)
                for word in words:
                    self.listbox.insert(tk.END, word)
            else:
                if self.listbox:
                    self.listbox.destroy()
                    self.listbox = None

    def selection(self, event):
        if self.listbox:
            self.var.set(self.listbox.get(tk.ACTIVE))
            self.listbox.destroy()
            self.listbox = None
            self.icursor(tk.END)

    def move_up(self, event):
        if self.listbox:
            if self.listbox.curselection() == (0,):
                return
            index = self.listbox.curselection()[0]
            self.listbox.selection_clear(first=index)
            index -= 1
            self.listbox.selection_set(first=index)
            self.listbox.activate(index)

    def move_down(self, event):
        if self.listbox:
            if self.listbox.curselection() == (self.listbox.size() - 1,):
                return
            index = self.listbox.curselection()[0]
            self.listbox.selection_clear(first=index)
            index += 1
            self.listbox.selection_set(first=index)
            self.listbox.activate(index)

    def comparison(self):
        return [w for w in self.autocompleteList if self.var.get().lower() in w.lower()]


def load_gb_file():
    filepath = filedialog.askopenfilename()
    if filepath:
        gb_file_path.set(filepath)
        update_annotations(filepath)

def update_annotations(gb_file_path):
    annotations = set()
    for record in SeqIO.parse(gb_file_path, "genbank"):
        for feature in record.features:
            if 'note' in feature.qualifiers:
                for note in feature.qualifiers['note']:
                    annotations.update(note.split())
    keyword_entry.update_autocompleteList(list(annotations))

def on_submit():
    selected_keyword = keyword_entry.get()
    selected_locus = locus_entry.get()
    selected_gb_file = gb_file_path.get()
    # Call your processing functions here using the selected parameters
    baseline_codon_counts, n_sequs = aggregate_codon_counts(selected_gb_file, selected_keyword)
    baseline_codon_usage = compute_codon_usage_from_counts(baseline_codon_counts)
    gene_codon_counts = get_codon_counts_for_gene_by_locus(selected_gb_file, selected_locus)
    gene_codon_usage = compute_codon_usage_from_counts(gene_codon_counts)
    compare_codon_usage(baseline_codon_usage, gene_codon_usage, n_sequs)
    validate_codon_usage_percentages(baseline_codon_usage)
    validate_codon_usage_percentages(gene_codon_usage)

root = tk.Tk()
root.title("Codon Usage Analysis GUI")

gb_file_path = tk.StringVar()

ttk.Label(root, text="Gene Locus:").pack(pady=5)
locus_entry = ttk.Entry(root)
locus_entry.pack()

ttk.Label(root, text="Select GB File:").pack(pady=5)
ttk.Entry(root, textvariable=gb_file_path, state='readonly').pack()
ttk.Button(root, text="Browse", command=load_gb_file).pack(pady=5)

ttk.Label(root, text="Enter Keyword:").pack(pady=5)
keyword_entry = AutocompleteEntry([], root, width=50)
keyword_entry.pack()

submit_btn = ttk.Button(root, text="Submit", command=on_submit)
submit_btn.pack(pady=10)

root.mainloop()
