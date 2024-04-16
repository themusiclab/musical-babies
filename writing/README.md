# How to use the additional files in the `/Writing` folder

As you write your manuscript, you have multiple citation style options to choose from. These styles are specified by .csl (Citation Style Language) files. While you might have to add a separate .csl file depending on the journal you're submitting to, the default in our repositories are either Nature or APA style. 

Notes on other customised files we keep here:

`apa-auto_cofirsts.csl` : This .csl file was put together by Courtney, and is an adjusted version of apa.csl that indicates co-first authors. For this file to work, update any 'citations.bib' file you're using such that for any citations you want to be shown as co-first authors, you add an additional 'note' field which contains the text "cofirst".

`bib_extractor.R` : Typically when you export a Zotero library, the .bib file will include all the information for every reference in that library. This is okay to do, but when we publish our data and analysis scripts online, it's cleaner to upload a "pruned" version of the .bib file that only includes information for the references being cited in the .rmd. There isn't an official way to do this, but this script takes your .rmd file and a .bib file, and spits out the pruned bibliography.