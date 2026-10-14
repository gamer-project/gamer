filename=spectral_tables.zip
link=https://use.yt/upload/66f39405

rm -r spectral_tables*
curl -L ${link} -o ${filename}
unzip ${filename}
