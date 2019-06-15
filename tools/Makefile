all: esmd

SRC_INC=$(subst :, -I, :$(VPATH))
include Depends
include $(TOOLROOT)/platforms/Macros.$(PLATFORM)

.PHONY: clean
clean:
	-rm -f $(OBJS) $(BLDDIR)/esmd
