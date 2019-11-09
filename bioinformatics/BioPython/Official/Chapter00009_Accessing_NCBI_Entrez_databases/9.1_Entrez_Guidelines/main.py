from Bio import Entrez

# ================================================================================
# @@ Before using Biopython to access the NCBI’s online resources (via Bio.Entrez or some of the other modules), please read the NCBI’s Entrez User Requirements.

# ================================================================================
# @@ If the NCBI finds you are abusing their systems, they can and will ban your access!

# ================================================================================
# @@ To paraphrase:

# ================================================================================
# @@ For any series of more than 100 requests, do this at weekends or outside USA peak times.

# ================================================================================
# @@ This is up to you to obey.

# ================================================================================
# @@ Use the https://eutils.ncbi.nlm.nih.gov address, not the standard NCBI Web address.

# ================================================================================
# @@ Biopython uses this web address.

# ================================================================================
# @@ If you are using a API key, you can make at most 10 queries per second, otherwise at most 3 queries per second.

# ================================================================================
# @@ This is automatically enforced by Biopython.

# ================================================================================
# Include api_key="MyAPIkey" in the argument list or set it as a module level variable:

# ================================================================================
# Entrez.api_key = "MyAPIkey"

# ================================================================================
# @@ Use the optional email parameter so the NCBI can contact you if there is a problem. 

# ================================================================================
# You can either explicitly set this as a parameter with each call to Entrez (e.g. include email="A.N.Other@example.com" in the argument list), or you can set a global email address:
# Entrez.email = "A.N.Other@example.com"

# ================================================================================
# @@ Bio.Entrez will then use this email address with each call to Entrez.

# ================================================================================
# @@ The example.com address is a reserved domain name specifically for documentation (RFC 2606).

# ================================================================================
# @@ Please DO NOT use a random email – it’s better not to give an email at all.

# ================================================================================
# @@ The email parameter has been mandatory since June 1, 2010.

# ================================================================================
# @@ In case of excessive usage, NCBI will attempt to contact a user at the e-mail address provided prior to blocking access to the E-utilities.

# ================================================================================
# @@ If you are using Biopython within some larger software suite, use the tool parameter to specify this.

# ================================================================================
# You can either explicitly set the tool name as a parameter with each call to Entrez (e.g. include tool="MyLocalScript" in the argument list), or you can set a global tool name:

# Entrez.tool = "MyLocalScript"

# ================================================================================
# The tool parameter will default to Biopython.

# ================================================================================
# For large queries, the NCBI also recommend using their session history feature (the WebEnv session cookie string, see Section 9.16). 

# ================================================================================
# This is only slightly more complicated.

# ================================================================================
# In conclusion, be sensible with your usage levels. 

# ================================================================================
# If you plan to download lots of data, consider other options. 

# ================================================================================
# For example, if you want easy access to all the human genes, consider fetching each chromosome by FTP as a GenBank file, and importing these into your own BioSQL database (see Section 20.5).

