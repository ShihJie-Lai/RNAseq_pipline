sub stringtie_RNA {
	my ($GFFGTF,$stringtieStrand,$stringtiePATH,$prep) = @_;
	$prep=~s/\/\w+.pl$//g;
	open(RSCRIPT,">stringTie_pip.sh") or die;
	print "\n$stringtiePATH/stringtie *.sorted.bam -p 16 -eB -x $stringtieStrand -o ./ballgown/\${aa}/\${aa}.gtf -G $GFFGTF -A ./ballgown/\${aa}/\${aa}_gene_abund.tab\n";
	print RSCRIPT"
	#!/bin/bash
	if [ ! -d ballgown ];then
		mkdir ballgown
	fi
	ls | grep .sorted.bam | grep -v '.bai' |while read id;
	do
		echo \$id
		aa=\$(echo \$id |sed  \"s/\.sorted\.bam//g\")
		if [ ! -d ./ballgown/\$aa ];then
			mkdir ./ballgown/\$aa
			$stringtiePATH/stringtie \${aa}.sorted.bam -p 16 -eB -x $stringtieStrand -o ./ballgown/\${aa}/\${aa}.gtf -G $GFFGTF -A ./ballgown/\${aa}/\${aa}_gene_abund.tab
		fi
		
	done
	python2.7 $prep/prepDE.py";

	close RSCRIPT;
	system("bash stringTie_pip.sh");
}
return 1;
