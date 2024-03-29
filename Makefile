#------------------------------------------------------------------------------#
# This makefile was generated by 'cbp2make' tool rev.147                       #
#------------------------------------------------------------------------------#


WORKDIR = `pwd`

CC = gcc
CXX = g++
AR = ar
LD = g++
WINDRES = windres

INC = 
CFLAGS = -Wall -fexceptions
RESINC = 
LIBDIR = 
LIB = libs/libjpeg.a libs/libpng.a libs/libz.a
LDFLAGS = 

INC_DEBUG = $(INC)
CFLAGS_DEBUG = $(CFLAGS) -g
RESINC_DEBUG = $(RESINC)
RCFLAGS_DEBUG = $(RCFLAGS)
LIBDIR_DEBUG = $(LIBDIR)
LIB_DEBUG = $(LIB)
LDFLAGS_DEBUG = $(LDFLAGS)
OBJDIR_DEBUG = obj/Debug
DEP_DEBUG = 
OUT_DEBUG = bin/Debug/lip-vireo

INC_RELEASE = $(INC)
CFLAGS_RELEASE = $(CFLAGS) -O3
RESINC_RELEASE = $(RESINC)
RCFLAGS_RELEASE = $(RCFLAGS)
LIBDIR_RELEASE = $(LIBDIR)
LIB_RELEASE = $(LIB)
LDFLAGS_RELEASE = $(LDFLAGS)
OBJDIR_RELEASE = obj/Release
DEP_RELEASE = 
OUT_RELEASE = bin/Release/lip-vireo

OBJ_DEBUG = $(OBJDIR_DEBUG)/iotool.o $(OBJDIR_DEBUG)/ioimage.o $(OBJDIR_DEBUG)/intimage.o $(OBJDIR_DEBUG)/image.o $(OBJDIR_DEBUG)/hesslap.o $(OBJDIR_DEBUG)/kernel.o $(OBJDIR_DEBUG)/hessian.o $(OBJDIR_DEBUG)/hessaff.o $(OBJDIR_DEBUG)/harris.o $(OBJDIR_DEBUG)/harlap.o $(OBJDIR_DEBUG)/haar.o $(OBJDIR_DEBUG)/filter.o $(OBJDIR_DEBUG)/pallete.o $(OBJDIR_DEBUG)/vstring.o $(OBJDIR_DEBUG)/vmath.o $(OBJDIR_DEBUG)/viewboard.o $(OBJDIR_DEBUG)/thinner.o $(OBJDIR_DEBUG)/scriptparser.o $(OBJDIR_DEBUG)/nondetector.o $(OBJDIR_DEBUG)/main.o $(OBJDIR_DEBUG)/log.o $(OBJDIR_DEBUG)/kpdrawer.o $(OBJDIR_DEBUG)/cimage.o $(OBJDIR_DEBUG)/desccm.o $(OBJDIR_DEBUG)/descasift.o $(OBJDIR_DEBUG)/dense.o $(OBJDIR_DEBUG)/cleaner.o  $(OBJDIR_DEBUG)/canny.o $(OBJDIR_DEBUG)/abstractimage.o $(OBJDIR_DEBUG)/abstractdetector.o $(OBJDIR_DEBUG)/abstractdescriptor.o $(OBJDIR_DEBUG)/descpview.o $(OBJDIR_DEBUG)/dsurf.o $(OBJDIR_DEBUG)/dog.o $(OBJDIR_DEBUG)/descsurf.o $(OBJDIR_DEBUG)/descspin.o $(OBJDIR_DEBUG)/descsift.o $(OBJDIR_DEBUG)/descrift.o $(OBJDIR_DEBUG)/descpcasift.o $(OBJDIR_DEBUG)/descljet.o $(OBJDIR_DEBUG)/descfind.o $(OBJDIR_DEBUG)/descfift.o $(OBJDIR_DEBUG)/descdelegator.o 

OBJ_RELEASE = $(OBJDIR_RELEASE)/iotool.o $(OBJDIR_RELEASE)/ioimage.o $(OBJDIR_RELEASE)/intimage.o $(OBJDIR_RELEASE)/image.o $(OBJDIR_RELEASE)/hesslap.o $(OBJDIR_RELEASE)/kernel.o $(OBJDIR_RELEASE)/hessian.o $(OBJDIR_RELEASE)/hessaff.o $(OBJDIR_RELEASE)/harris.o $(OBJDIR_RELEASE)/harlap.o $(OBJDIR_RELEASE)/haar.o $(OBJDIR_RELEASE)/filter.o $(OBJDIR_RELEASE)/pallete.o $(OBJDIR_RELEASE)/vstring.o $(OBJDIR_RELEASE)/vmath.o $(OBJDIR_RELEASE)/viewboard.o $(OBJDIR_RELEASE)/thinner.o $(OBJDIR_RELEASE)/scriptparser.o $(OBJDIR_RELEASE)/nondetector.o $(OBJDIR_RELEASE)/main.o $(OBJDIR_RELEASE)/log.o $(OBJDIR_RELEASE)/kpdrawer.o $(OBJDIR_RELEASE)/cimage.o $(OBJDIR_RELEASE)/desccm.o $(OBJDIR_RELEASE)/descasift.o $(OBJDIR_RELEASE)/dense.o $(OBJDIR_RELEASE)/cleaner.o $(OBJDIR_RELEASE)/canny.o $(OBJDIR_RELEASE)/abstractimage.o $(OBJDIR_RELEASE)/abstractdetector.o $(OBJDIR_RELEASE)/abstractdescriptor.o $(OBJDIR_RELEASE)/descpview.o $(OBJDIR_RELEASE)/dsurf.o $(OBJDIR_RELEASE)/dog.o $(OBJDIR_RELEASE)/descsurf.o $(OBJDIR_RELEASE)/descspin.o $(OBJDIR_RELEASE)/descsift.o $(OBJDIR_RELEASE)/descrift.o  $(OBJDIR_RELEASE)/descpcasift.o $(OBJDIR_RELEASE)/descljet.o $(OBJDIR_RELEASE)/descfind.o $(OBJDIR_RELEASE)/descfift.o $(OBJDIR_RELEASE)/descdelegator.o

all: debug release

clean: clean_debug clean_release

before_debug: 
	test -d bin/Debug || mkdir -p bin/Debug
	test -d $(OBJDIR_DEBUG) || mkdir -p $(OBJDIR_DEBUG)

after_debug: 

debug: before_debug out_debug after_debug

out_debug: before_debug $(OBJ_DEBUG) $(DEP_DEBUG)
	$(LD) $(LIBDIR_DEBUG) -o $(OUT_DEBUG) $(OBJ_DEBUG)  $(LDFLAGS_DEBUG) $(LIB_DEBUG)

$(OBJDIR_DEBUG)/iotool.o: iotool.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c iotool.cpp -o $(OBJDIR_DEBUG)/iotool.o

$(OBJDIR_DEBUG)/ioimage.o: ioimage.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c ioimage.cpp -o $(OBJDIR_DEBUG)/ioimage.o

$(OBJDIR_DEBUG)/intimage.o: intimage.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c intimage.cpp -o $(OBJDIR_DEBUG)/intimage.o

$(OBJDIR_DEBUG)/image.o: image.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c image.cpp -o $(OBJDIR_DEBUG)/image.o

$(OBJDIR_DEBUG)/hesslap.o: hesslap.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c hesslap.cpp -o $(OBJDIR_DEBUG)/hesslap.o

$(OBJDIR_DEBUG)/kernel.o: kernel.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c kernel.cpp -o $(OBJDIR_DEBUG)/kernel.o

$(OBJDIR_DEBUG)/hessian.o: hessian.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c hessian.cpp -o $(OBJDIR_DEBUG)/hessian.o

$(OBJDIR_DEBUG)/hessaff.o: hessaff.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c hessaff.cpp -o $(OBJDIR_DEBUG)/hessaff.o

$(OBJDIR_DEBUG)/harris.o: harris.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c harris.cpp -o $(OBJDIR_DEBUG)/harris.o

$(OBJDIR_DEBUG)/harlap.o: harlap.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c harlap.cpp -o $(OBJDIR_DEBUG)/harlap.o

$(OBJDIR_DEBUG)/haar.o: haar.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c haar.cpp -o $(OBJDIR_DEBUG)/haar.o

$(OBJDIR_DEBUG)/filter.o: filter.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c filter.cpp -o $(OBJDIR_DEBUG)/filter.o

$(OBJDIR_DEBUG)/pallete.o: pallete.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c pallete.cpp -o $(OBJDIR_DEBUG)/pallete.o

$(OBJDIR_DEBUG)/vstring.o: vstring.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c vstring.cpp -o $(OBJDIR_DEBUG)/vstring.o

$(OBJDIR_DEBUG)/vmath.o: vmath.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c vmath.cpp -o $(OBJDIR_DEBUG)/vmath.o

$(OBJDIR_DEBUG)/viewboard.o: viewboard.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c viewboard.cpp -o $(OBJDIR_DEBUG)/viewboard.o

$(OBJDIR_DEBUG)/thinner.o: thinner.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c thinner.cpp -o $(OBJDIR_DEBUG)/thinner.o

$(OBJDIR_DEBUG)/scriptparser.o: scriptparser.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c scriptparser.cpp -o $(OBJDIR_DEBUG)/scriptparser.o

$(OBJDIR_DEBUG)/nondetector.o: nondetector.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c nondetector.cpp -o $(OBJDIR_DEBUG)/nondetector.o

$(OBJDIR_DEBUG)/main.o: main.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c main.cpp -o $(OBJDIR_DEBUG)/main.o

$(OBJDIR_DEBUG)/log.o: log.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c log.cpp -o $(OBJDIR_DEBUG)/log.o

$(OBJDIR_DEBUG)/kpdrawer.o: kpdrawer.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c kpdrawer.cpp -o $(OBJDIR_DEBUG)/kpdrawer.o

$(OBJDIR_DEBUG)/cimage.o: cimage.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c cimage.cpp -o $(OBJDIR_DEBUG)/cimage.o

$(OBJDIR_DEBUG)/desccm.o: desccm.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c desccm.cpp -o $(OBJDIR_DEBUG)/desccm.o

$(OBJDIR_DEBUG)/descasift.o: descasift.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c descasift.cpp -o $(OBJDIR_DEBUG)/descasift.o

$(OBJDIR_DEBUG)/dense.o: dense.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c dense.cpp -o $(OBJDIR_DEBUG)/dense.o

$(OBJDIR_DEBUG)/cleaner.o: cleaner.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c cleaner.cpp -o $(OBJDIR_DEBUG)/cleaner.o

$(OBJDIR_DEBUG)/canny.o: canny.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c canny.cpp -o $(OBJDIR_DEBUG)/canny.o

$(OBJDIR_DEBUG)/abstractimage.o: abstractimage.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c abstractimage.cpp -o $(OBJDIR_DEBUG)/abstractimage.o

$(OBJDIR_DEBUG)/abstractdetector.o: abstractdetector.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c abstractdetector.cpp -o $(OBJDIR_DEBUG)/abstractdetector.o

$(OBJDIR_DEBUG)/abstractdescriptor.o: abstractdescriptor.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c abstractdescriptor.cpp -o $(OBJDIR_DEBUG)/abstractdescriptor.o

$(OBJDIR_DEBUG)/descpview.o: descpview.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c descpview.cpp -o $(OBJDIR_DEBUG)/descpview.o

$(OBJDIR_DEBUG)/dsurf.o: dsurf.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c dsurf.cpp -o $(OBJDIR_DEBUG)/dsurf.o

$(OBJDIR_DEBUG)/dog.o: dog.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c dog.cpp -o $(OBJDIR_DEBUG)/dog.o

$(OBJDIR_DEBUG)/descsurf.o: descsurf.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c descsurf.cpp -o $(OBJDIR_DEBUG)/descsurf.o

$(OBJDIR_DEBUG)/descspin.o: descspin.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c descspin.cpp -o $(OBJDIR_DEBUG)/descspin.o

$(OBJDIR_DEBUG)/descsift.o: descsift.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c descsift.cpp -o $(OBJDIR_DEBUG)/descsift.o

$(OBJDIR_DEBUG)/descrift.o: descrift.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c descrift.cpp -o $(OBJDIR_DEBUG)/descrift.o

$(OBJDIR_DEBUG)/descpcasift.o: descpcasift.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c descpcasift.cpp -o $(OBJDIR_DEBUG)/descpcasift.o

$(OBJDIR_DEBUG)/descljet.o: descljet.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c descljet.cpp -o $(OBJDIR_DEBUG)/descljet.o

$(OBJDIR_DEBUG)/descfind.o: descfind.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c descfind.cpp -o $(OBJDIR_DEBUG)/descfind.o

$(OBJDIR_DEBUG)/descfift.o: descfift.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c descfift.cpp -o $(OBJDIR_DEBUG)/descfift.o

$(OBJDIR_DEBUG)/descdelegator.o: descdelegator.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c descdelegator.cpp -o $(OBJDIR_DEBUG)/descdelegator.o

$(OBJDIR_DEBUG)/desccontour.o: desccontour.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c desccontour.cpp -o $(OBJDIR_DEBUG)/desccontour.o

clean_debug: 
	rm -f $(OBJ_DEBUG) $(OUT_DEBUG)
	rm -rf bin/Debug
	rm -rf $(OBJDIR_DEBUG)

before_release: 
	test -d bin/Release || mkdir -p bin/Release
	test -d $(OBJDIR_RELEASE) || mkdir -p $(OBJDIR_RELEASE)

after_release: 

release: before_release out_release after_release

out_release: before_release $(OBJ_RELEASE) $(DEP_RELEASE)
	$(LD) $(LIBDIR_RELEASE) -o $(OUT_RELEASE) $(OBJ_RELEASE)  $(LDFLAGS_RELEASE) $(LIB_RELEASE)

$(OBJDIR_RELEASE)/iotool.o: iotool.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c iotool.cpp -o $(OBJDIR_RELEASE)/iotool.o

$(OBJDIR_RELEASE)/ioimage.o: ioimage.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c ioimage.cpp -o $(OBJDIR_RELEASE)/ioimage.o

$(OBJDIR_RELEASE)/intimage.o: intimage.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c intimage.cpp -o $(OBJDIR_RELEASE)/intimage.o

$(OBJDIR_RELEASE)/image.o: image.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c image.cpp -o $(OBJDIR_RELEASE)/image.o

$(OBJDIR_RELEASE)/hesslap.o: hesslap.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c hesslap.cpp -o $(OBJDIR_RELEASE)/hesslap.o

$(OBJDIR_RELEASE)/kernel.o: kernel.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c kernel.cpp -o $(OBJDIR_RELEASE)/kernel.o

$(OBJDIR_RELEASE)/hessian.o: hessian.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c hessian.cpp -o $(OBJDIR_RELEASE)/hessian.o

$(OBJDIR_RELEASE)/hessaff.o: hessaff.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c hessaff.cpp -o $(OBJDIR_RELEASE)/hessaff.o

$(OBJDIR_RELEASE)/harris.o: harris.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c harris.cpp -o $(OBJDIR_RELEASE)/harris.o

$(OBJDIR_RELEASE)/harlap.o: harlap.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c harlap.cpp -o $(OBJDIR_RELEASE)/harlap.o

$(OBJDIR_RELEASE)/haar.o: haar.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c haar.cpp -o $(OBJDIR_RELEASE)/haar.o

$(OBJDIR_RELEASE)/filter.o: filter.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c filter.cpp -o $(OBJDIR_RELEASE)/filter.o

$(OBJDIR_RELEASE)/pallete.o: pallete.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c pallete.cpp -o $(OBJDIR_RELEASE)/pallete.o

$(OBJDIR_RELEASE)/vstring.o: vstring.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c vstring.cpp -o $(OBJDIR_RELEASE)/vstring.o

$(OBJDIR_RELEASE)/vmath.o: vmath.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c vmath.cpp -o $(OBJDIR_RELEASE)/vmath.o

$(OBJDIR_RELEASE)/viewboard.o: viewboard.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c viewboard.cpp -o $(OBJDIR_RELEASE)/viewboard.o

$(OBJDIR_RELEASE)/thinner.o: thinner.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c thinner.cpp -o $(OBJDIR_RELEASE)/thinner.o

$(OBJDIR_RELEASE)/scriptparser.o: scriptparser.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c scriptparser.cpp -o $(OBJDIR_RELEASE)/scriptparser.o

$(OBJDIR_RELEASE)/nondetector.o: nondetector.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c nondetector.cpp -o $(OBJDIR_RELEASE)/nondetector.o

$(OBJDIR_RELEASE)/main.o: main.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c main.cpp -o $(OBJDIR_RELEASE)/main.o

$(OBJDIR_RELEASE)/log.o: log.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c log.cpp -o $(OBJDIR_RELEASE)/log.o

$(OBJDIR_RELEASE)/kpdrawer.o: kpdrawer.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c kpdrawer.cpp -o $(OBJDIR_RELEASE)/kpdrawer.o

$(OBJDIR_RELEASE)/cimage.o: cimage.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c cimage.cpp -o $(OBJDIR_RELEASE)/cimage.o

$(OBJDIR_RELEASE)/desccm.o: desccm.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c desccm.cpp -o $(OBJDIR_RELEASE)/desccm.o

$(OBJDIR_RELEASE)/descasift.o: descasift.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c descasift.cpp -o $(OBJDIR_RELEASE)/descasift.o

$(OBJDIR_RELEASE)/dense.o: dense.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c dense.cpp -o $(OBJDIR_RELEASE)/dense.o

$(OBJDIR_RELEASE)/cleaner.o: cleaner.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c cleaner.cpp -o $(OBJDIR_RELEASE)/cleaner.o

$(OBJDIR_RELEASE)/canny.o: canny.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c canny.cpp -o $(OBJDIR_RELEASE)/canny.o

$(OBJDIR_RELEASE)/abstractimage.o: abstractimage.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c abstractimage.cpp -o $(OBJDIR_RELEASE)/abstractimage.o

$(OBJDIR_RELEASE)/abstractdetector.o: abstractdetector.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c abstractdetector.cpp -o $(OBJDIR_RELEASE)/abstractdetector.o

$(OBJDIR_RELEASE)/abstractdescriptor.o: abstractdescriptor.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c abstractdescriptor.cpp -o $(OBJDIR_RELEASE)/abstractdescriptor.o

$(OBJDIR_RELEASE)/descpview.o: descpview.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c descpview.cpp -o $(OBJDIR_RELEASE)/descpview.o

$(OBJDIR_RELEASE)/dsurf.o: dsurf.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c dsurf.cpp -o $(OBJDIR_RELEASE)/dsurf.o

$(OBJDIR_RELEASE)/dog.o: dog.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c dog.cpp -o $(OBJDIR_RELEASE)/dog.o

$(OBJDIR_RELEASE)/descsurf.o: descsurf.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c descsurf.cpp -o $(OBJDIR_RELEASE)/descsurf.o

$(OBJDIR_RELEASE)/descspin.o: descspin.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c descspin.cpp -o $(OBJDIR_RELEASE)/descspin.o

$(OBJDIR_RELEASE)/descsift.o: descsift.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c descsift.cpp -o $(OBJDIR_RELEASE)/descsift.o

$(OBJDIR_RELEASE)/descrift.o: descrift.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c descrift.cpp -o $(OBJDIR_RELEASE)/descrift.o

$(OBJDIR_RELEASE)/descpcasift.o: descpcasift.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c descpcasift.cpp -o $(OBJDIR_RELEASE)/descpcasift.o

$(OBJDIR_RELEASE)/descljet.o: descljet.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c descljet.cpp -o $(OBJDIR_RELEASE)/descljet.o

$(OBJDIR_RELEASE)/descfind.o: descfind.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c descfind.cpp -o $(OBJDIR_RELEASE)/descfind.o

$(OBJDIR_RELEASE)/descfift.o: descfift.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c descfift.cpp -o $(OBJDIR_RELEASE)/descfift.o

$(OBJDIR_RELEASE)/descdelegator.o: descdelegator.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c descdelegator.cpp -o $(OBJDIR_RELEASE)/descdelegator.o

$(OBJDIR_RELEASE)/desccontour.o: desccontour.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c desccontour.cpp -o $(OBJDIR_RELEASE)/desccontour.o

clean_release: 
	rm -f $(OBJ_RELEASE) $(OUT_RELEASE)
	rm -rf bin/Release
	rm -rf $(OBJDIR_RELEASE)

.PHONY: before_debug after_debug clean_debug before_release after_release clean_release

