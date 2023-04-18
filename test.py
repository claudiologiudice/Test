from Bio import Entrez
from Bio import SeqIO


# Imposta i parametri di ricerca
query = "(misc_feature OR cds OR rrna OR trna) AND (\"rna editing\" OR \"rna edited\") AND (\"plastid\"[Filter] OR \"mitochondrion\"[Filter]) AND \"refseq\"[Filter] AND \"plants\"[Organism]"

# Imposta le credenziali di accesso a Entrez
Entrez.email = ""
Entrez.api_key = ""


# Esegui la ricerca su GenBank
handle = Entrez.esearch(db="nucleotide", term=query, retmax=1000)
record_ids = Entrez.read(handle)["IdList"]

# Scarica i record da GenBank
records = []
for record_id in record_ids:
	handle = Entrez.efetch(db="nucleotide", id=record_id, rettype="gb", retmode="text")
	record = SeqIO.read(handle, "genbank")
	records.append(record)
	# Applica lo script a ciascun record
	for record in records:
	   # Crea un dizionario vuoto per contenere le informazioni delle CDS
		cds_dict = {}

		# Crea un dizionario vuoto per contenere le informazioni degli rRNA
		rrna_dict = {}

		# Crea un dizionario vuoto per contenere le informazioni dei tRNA
		trna_dict = {}

		# Itera attraverso tutte le feature del record e seleziona solo le CDS
		for feature in record.features:
			if feature.type == "CDS":
				cds_dict.setdefault((feature.location.start, feature.location.end), []) 

		# Itera attraverso tutte le feature del record e seleziona solo gli rRNA
		for feature in record.features:
			if feature.type == "rRNA":
				rrna_dict.setdefault((feature.location.start, feature.location.end), []) 

		# Itera attraverso tutte le feature del record e seleziona solo i tRNA
		for feature in record.features:
			if feature.type == "tRNA":
				trna_dict.setdefault((feature.location.start, feature.location.end), []) 

		# Itera nuovamente attraverso tutte le feature del record e seleziona solo le misc_feature che contengono "rna editing" o "rna e" nelle 
note
		for feature in record.features:
			if feature.type == "misc_feature":
				if "rna editing" in feature.qualifiers.get("note", [""])[0].lower() or "rna e" in feature.qualifiers.get("note", [""])[0].lower():
					for cds in cds_dict:
						if cds[0] <= feature.location.start <= cds[1] or cds[0] <= feature.location.end <= cds[1]:
							cds_dict[cds].append(feature)
					for rrna in rrna_dict:
						if rrna[0] <= feature.location.start <= rrna[1] or rrna[0] <= feature.location.end <= rrna[1]:
							rrna_dict[rrna].append(feature)
					for trna in trna_dict:
						if trna[0] <= feature.location.start <= trna[1] or trna[0] <= feature.location.end <= trna[1]:
							trna_dict[trna].append(feature)

		# Stampa i dizionari
		if any(len(value) > 0 for value in cds_dict.values()):
				print(record.id)
				print("CDS:")
				print(cds_dict)
		if any(len(value) > 0 for value in rrna_dict.values()):
				print(record.id)
				print("rRNA:")
				print(rrna_dict)
		if any(len(value) > 0 for value in trna_dict.values()):
				print(record.id)
				print("tRNA:")
				print(trna_dict)
