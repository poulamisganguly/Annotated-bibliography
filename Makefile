PAPERS = $(wildcard papers/*.md)
HTMLS  = $(subst papers/,html/,$(patsubst %.md,%.html,$(PAPERS)))

default: $(HTMLS)

html/%.html: papers/%.md header.html footer.html
	cat header.html > $@
	pandoc --mathjax --toc --from markdown --to html --ascii $< >> $@
	cat footer.html >> $@
