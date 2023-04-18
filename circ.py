import requests
from Bio import Entrez, SeqIO
from jinja2 import Template

# specifica le tue credenziali di accesso a Entrez
Entrez.email = "claudio.logiudice@uniba.it"

# cerca il record con l'ID desiderato
handle = Entrez.efetch(db="nucleotide", id="AB501189", rettype="gb", retmode="text")

# Carica il file FASTA del genoma
genome = SeqIO.read(handle, "gb")

# Estrae i nomi dei geni e le loro posizioni
genes = []
for feature in genome.features:
    if feature.type == "gene":
        genes.append({
            "name": feature.qualifiers.get("gene", [""])[0],
            "start": feature.location.start.position,
            "end": feature.location.end.position
        })

# Ordina i geni per posizione di inizio
genes = sorted(genes, key=lambda gene: gene["start"])
