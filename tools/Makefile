all: esmd

SRC_INC=$(subst :, -I, :$(VPATH))
include Depends
include $(CASEROOT)/Macros

.PHONY: clean
clean:
	-rm -f $(OBJS) $(BLDDIR)/esmd
