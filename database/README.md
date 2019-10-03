# Database for Spyctra

Prototyping the Spyctra database

Spyctra will be shipped with a single data file (besides the Vega spectrum) which
will contain the web location of the database

The database will be organized as follow (re-using a bit the synphot CDBS structure) 


```
database
│   index.html
│   database_contents.yaml?    
│
└───templates
│   │
│   └─── library1
│   │   │  contents.yaml?
│   │   │   spectral_template1.fits
│   │   │   spectral_template2.fits
│   │   │   ...
│   │ 
│   └─── library2
│       │   contents.yaml?
│       │   another_spectral_template1.ascii
│       │    ...
│   
└───extinction
│   │   contents.yml?
│   │   extinction_curve1.fits
│   │   extinction_curve2.fits
│   │   ...
│
└───filters?          
```   

Questions marks indicate that no decision has been made yet. 

* At the root directory, there should be a file describing the contents of 
the database. Every time that data is added this file should be updated.

* Ideally the format of the data file should allow to be displayed in a browser for 
a quick look.  Not sure if yaml fits the bill (hence the question marks)

* Similarly each data directory should contain a file describing the contents, 
like list of templates, a summary and references.

* Data files should be readable by synphot

* **Content files should be able to be parsed by Spyctra and display the contents 
to the user and/or pass them to the API**

Following the same standard other datafiles can be added like filters, etc. 

