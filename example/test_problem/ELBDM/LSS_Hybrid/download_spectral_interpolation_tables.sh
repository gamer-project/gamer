filename=spectral_tables.zip
link=http://use.yt/upload/66f39405

rm -r spectral_tables*
curl ${link} -o ${filename}
unzip ${filename}