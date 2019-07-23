all: esmd

SRC_INC=$(subst :, -I, :$(SRCPATH))
include Depends
include $(CASEROOT)/Macros

.PHONY: clean
clean:
	-rm -f $(OBJS) $(BLDDIR)/esmd
