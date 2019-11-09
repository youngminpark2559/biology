# Entrez (https://www.ncbi.nlm.nih.gov/Web/Search/entrezfs.html) is a data retrieval system that provides users access to NCBI’s databases such as PubMed, GenBank, GEO, and many others.

# users --> Entrez --> NCBI’s databases (PubMed, GenBank, GEO, and many others)

# ================================================================================
# @@ You can access Entrez from a web browser to manually enter queries, or you can use Biopython’s Bio.Entrez module for programmatic access to Entrez.

# ================================================================================
# The latter allows you for example to search PubMed or download GenBank records from within a Python script.

# PubMed or download GenBank records = Bio.Entrez(your request)

# ================================================================================
# The Bio.Entrez module makes use of the Entrez Programming Utilities (also known as EUtils), consisting of eight tools that are described in detail on NCBI’s page at https://www.ncbi.nlm.nih.gov/books/NBK25501/.

# Bio.Entrez uses EUtils (Entrez Programming Utilities)

# EUtils = 8 tools

# ================================================================================
# @@ Each of these tools corresponds to one Python function in the Bio.Entrez module, as described in the sections below.

# ================================================================================
# @@ This module makes sure that the correct URL is used for the queries, and that not more than one request is made every three seconds, as required by NCBI.

# ================================================================================
# The output returned by the Entrez Programming Utilities is typically in XML format.

# ================================================================================
# To parse such output, you have several options:

# ================================================================================
# Use Bio.Entrez’s parser to parse the XML output into a Python object;

# Python object = Bio.Entrez’s parser(XML output)

# ================================================================================
# Use the DOM (Document Object Model) parser in Python’s standard library;

# parsed output = Document Object Model parser(XML output)

# ================================================================================
# Use the SAX (Simple API for XML) parser in Python’s standard library;

# parsed output = Simple API for XML(XML output)

# ================================================================================
# @@ Read the XML output as raw text, and parse it by string searching and manipulation.

# ================================================================================
# @@ For the DOM and SAX parsers, see the Python documentation.

# ================================================================================
# @@ The parser in Bio.Entrez is discussed below.

# ================================================================================
# @@ NCBI uses DTD (Document Type Definition) files to describe the structure of the information contained in XML files.

# ================================================================================
# @@ Most of the DTD files used by NCBI are included in the Biopython distribution.

# ================================================================================
# @@ The Bio.Entrez parser makes use of the DTD files when parsing an XML file returned by NCBI Entrez.

# ================================================================================
# @@ Occasionally, you may find that the DTD file associated with a specific XML file is missing in the Biopython distribution.

# ================================================================================
# @@ In particular, this may happen when NCBI updates its DTD files.

# ================================================================================
# @@ If this happens, Entrez.read will show a warning message with the name and URL of the missing DTD file.

# ================================================================================
# @@ The parser will proceed to access the missing DTD file through the internet, allowing the parsing of the XML file to continue.

# ================================================================================
# However, the parser is much faster if the DTD file is available locally.

# ================================================================================
# For this purpose, please download the DTD file from the URL in the warning message and place it in the directory ...site-packages/Bio/Entrez/DTDs, containing the other DTD files.

# ================================================================================
# If you don’t have write access to this directory, you can also place the DTD file in ~/.biopython/Bio/Entrez/DTDs, where ~ represents your home directory.

# ================================================================================
# Since this directory is read before the directory ...site-packages/Bio/Entrez/DTDs, you can also put newer versions of DTD files there if the ones in ...site-packages/Bio/Entrez/DTDs become outdated.

# ================================================================================
# Alternatively, if you installed Biopython from source, you can add the DTD file to the source code’s Bio/Entrez/DTDs directory, and reinstall Biopython.

# ================================================================================
# This will install the new DTD file in the correct location together with the other DTD files.

# ================================================================================
# The Entrez Programming Utilities can also generate output in other formats, such as the Fasta or GenBank file formats for sequence databases, or the MedLine format for the literature database, discussed in Section 9.13.
