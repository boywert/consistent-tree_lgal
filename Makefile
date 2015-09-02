SUBDIRS = src lgal_convert

.PHONY: subdirs $(SUBDIRS)

subdirs: $(SUBDIRS) 

all: subdirs

debug:
	-for d in $(SUBDIRS); do (cd $$d; $(MAKE) debug ); done

no-periodic:
	-for d in $(SUBDIRS); do (cd $$d; $(MAKE) no-periodic ); done

prof:
	-for d in $(SUBDIRS); do (cd $$d; $(MAKE) prof ); done

clean:
	-for d in $(SUBDIRS); do (cd $$d; $(MAKE) clean ); done
	-rm -f *~

.REMAKE:

dist: .REMAKE
	@perl -ne 'print "$$1\n" if (/VERSION\s*\"([^\"]+)/)' src/version.h > VERSION
	cd ../ ; tar -czvf ctrees.tar.gz consistent_trees/Makefile consistent_trees/*.pl consistent_trees/*.cfg consistent_trees/*/*.[ch] consistent_trees/*/Makefile consistent_trees/*/*.pm consistent_trees/src/*/*.pm consistent_trees/README.pdf  consistent_trees/README.md  consistent_trees/LICENSE consistent_trees/VERSION consistent_trees/CHANGELOG; mv ctrees.tar.gz consistent_trees

versiondist:
	rm -rf dist
	mkdir dist
	cd dist; tar xzf ../ctrees.tar.gz ; perl -ne '/\#define.*VERSION\D*([\d\.]+)/ && print $$1' consistent_trees/src/version.h > NUMBER ; mv consistent_trees consistent_trees-`cat NUMBER`; tar czf ctrees-`cat NUMBER`.tar.gz consistent_trees-`cat NUMBER`

$(SUBDIRS):
	$(MAKE) -C $@
	make no-periodic
