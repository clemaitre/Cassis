# Installation

### Lastz

[lastz](https://github.com/lastz/lastz)

Sur genocluster :

```
ll /local/env/envlastz-1.0.4.sh 
-rw-r--r-- 1 root root 80  4 févr.  2021 /local/env/envlastz-1.0.4.sh
```

prebuilt for MacOs : [here](http://hgdownload.soe.ucsc.edu/admin/exe/macOSX.x86_64/)
installé dans `/Users/clemaitr/Bin/lastz-1.04.00`

```
cd lastz/
ln -s /Users/clemaitr/Bin/lastz-1.04.00 lastz
./lastz
```

A l'air de marcher.


lastz and lastz_D
The two executables are basically the same program; the only difference is that lastz uses integer scores, while lastz_D uses floating-point scores. 

### Ghostscript

Teste de la commande bitmap :

```
R
bitmap("test.png",type="png16m", height=8, width=7, res=100)
plot(rnorm(10),rnorm(10))
dev.off()
```

ca marche !


## Sur genocluster

```
. /local/env/envlastz-1.0.4.sh 
. /local/env/envr-3.6.2.sh
. /local/env/envghostscript-9.53.3.sh

whereis lastz
cd Cassis/lastz
ln -s /softs/local/miniconda3/envs/lastz-1.0.4/bin/lastz .

cd ../test_papillons
mkdir res_cassis
perl ../cassis.pl BST1_LAU8.scaffold_4.asynt.nucmer.blocks.tab B BST1 LAU8 res_cassis
```


# A modifier

- paramètre `--lastzlevel`ne marche plus car l'option à changer dans lastz : 
	- modifier `core/alignSequences.pl` l. 152 
	
	```
	$command .= " Q=${matrix} ";   -> $command .= " --scores=${matrix} ";
	```
	- test : n'a pas marché... alignements idem...
	
	
