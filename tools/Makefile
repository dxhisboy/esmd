all: $(BLDDIR)/esmd

SRC_INC=$(subst :, -I, :$(VPATH))
include Depends
include $(TOOLROOT)/Macros/Macros.$(PLATFORM)

$(BLDDIR)/esmd: $(HISTDIR)/esmd.$(LID)
	cp $< $@
$(HISTDIR)/esmd.$(LID): esmd
	cp $< $@
