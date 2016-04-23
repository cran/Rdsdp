
default: 
	make all

D1       =  src
D2       =  lib

include make.include

dsdplibrary: chkroot info
	cd ${D1}; make all

clean:
	cd ${D1}; make clean
	cd ${D2}; make clean
	${RM} *~ */*~ */*/*~

htmlzip:
	zip -r DSDP5-api-html.zip dox
	${RM} -R dox

oshared: 
	-@${RM} tmp; \
	mkdir tmp; \
	cd tmp; \
	echo "building ${DSDPLIBSO}"; \
	${AR} x ${DSDPLIB} ;\
	${SH_LD} ${DSDPLIBSO} *.o -o ${DSDPLIBSO}; \
	cd ../ ; \
	${RM} tmp
