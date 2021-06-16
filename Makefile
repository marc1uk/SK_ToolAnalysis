include $(SKOFL_ROOT)/config.gmk  # pulls in libskroot.so as well

# user library paths
LDFLAGS += -L${HOME}/skrootlibs -L${HOME}/stllibs -L${HOME}/relic_sk4_ana/relic_work_dir/data_reduc/third/lib
# user libraries
LOCAL_LIBS = -lRootStl -lthirdredvars

# lowe libraries - some of these may not be required in this list
LDLIBS += -lbonsai_3.3 -lsklowe_7.0 -lwtlib_5.1 -lsollib_4.0 -lskrd -lsklib -lskroot -liolib -llibrary
LDLIBS += `cernlib graflib grafX11 packlib mathlib kernlib lapack3 blas`

LD_RUN_PATH=$(SKOFL_LIBDIR):$(A_LIBDIR)

# C++ compiler flags - XXX config.gmk sets this already, so APPEND ONLY XXX
CXXFLAGS    += -g -O3 -std=c++11 -fdiagnostics-color=always -Wno-reorder -Wno-sign-compare -Wno-unused-variable -Wno-unused-but-set-variable -Werror=array-bounds
#-D_GLIBCXX_DEBUG  << g++ debug mode. immediate segfault...

# flags required for gprof profiling
#CXXFLAGS    += -g -pg -ggdb3

# ToolDAQFramework debug mode: disable the try{}-catch{} around all Tool methods.
# Combine with -lSegFault to cause exceptions to invoke a segfault, printing a backtrace.
#CXXFLAGS     += -DDEBUG
#LDFLAGS      += -lSegFault

ToolDAQPath=ToolDAQ
ZMQLib= -L $(ToolDAQPath)/zeromq-4.0.7/lib -lzmq 
ZMQInclude= -I $(ToolDAQPath)/zeromq-4.0.7/include/ 

BoostLib= -L $(ToolDAQPath)/boost_1_66_0/install/lib -lboost_date_time -lboost_serialization -lboost_iostreams -lboost_regex
BoostInclude= -I $(ToolDAQPath)/boost_1_66_0/install/include

RootInclude= `root-config --cflags`
RootLib= `root-config --libs --evelibs --glibs`
# --glibs
# --evelibs for TParticlePDG

PythonInclude= `python3-config --cflags`
PythonLib = `python3-config --libs`

DataModelInclude = $(RootInclude)
DataModelLib = $(RootLib)

MyToolsInclude = $(PythonInclude)
MyToolsLib = $(LDFLAGS) $(LDLIBS) $(PythonLib)

all: lib/libStore.so lib/libLogging.so lib/libDataModel.so include/Tool.h lib/libMyTools.so lib/libServiceDiscovery.so lib/libToolChain.so main RemoteControl  NodeDaemon

main: src/main.cpp | lib/libMyTools.so lib/libStore.so lib/libLogging.so lib/libToolChain.so lib/libDataModel.so lib/libServiceDiscovery.so lib/liblowfit_sk4_stripped.so
	@echo -e "\n*************** Making " $@ "****************"
	g++ $(CXXFLAGS) -g -L lib -llowfit_sk4_stripped -I include $(DataModelInclude) $(BoostInclude) $(ZMQInclude) $(MyToolsInclude) src/main.cpp -o $@ $(BoostLib) $(DataModelLib) $(MyToolsLib) $(ZMQLib) -L lib -lStore -lMyTools -lToolChain -lDataModel -lLogging -lServiceDiscovery -lpthread

lib/libStore.so: $(ToolDAQPath)/ToolDAQFramework/src/Store/*
	cd $(ToolDAQPath)/ToolDAQFramework && make lib/libStore.so
	@echo -e "\n*************** Copying " $@ "****************"
	cp $(ToolDAQPath)/ToolDAQFramework/src/Store/*.h include/
	cp $(ToolDAQPath)/ToolDAQFramework/lib/libStore.so lib/
	#g++ -g -fPIC -shared  -I include $(ToolDAQPath)/ToolDAQFramework/src/Store/*.cpp -o lib/libStore.so $(BoostLib) $(BoostInclude)


include/Tool.h:  $(ToolDAQPath)/ToolDAQFramework/src/Tool/Tool.h
	@echo -e "\n*************** Copying " $@ "****************"
	cp $(ToolDAQPath)/ToolDAQFramework/src/Tool/Tool.h include/
	cp UserTools/*.h include/
	cp UserTools/*/*.h include/
	cp DataModel/*.h include/


lib/libToolChain.so: $(ToolDAQPath)/ToolDAQFramework/src/ToolChain/* | lib/libLogging.so lib/libStore.so lib/libMyTools.so lib/libServiceDiscovery.so lib/libLogging.so lib/libDataModel.so
	@echo -e "\n*************** Making " $@ "****************"
	cp $(ToolDAQPath)/ToolDAQFramework/UserTools/Factory/*.h include/
	cp $(ToolDAQPath)/ToolDAQFramework/src/ToolChain/*.h include/
	g++ $(CXXFLAGS) -g -fPIC -shared $(ToolDAQPath)/ToolDAQFramework/src/ToolChain/ToolChain.cpp -I include -lpthread -L lib -lStore -lDataModel -lServiceDiscovery -lLogging -lMyTools -o lib/libToolChain.so $(DataModelInclude) $(DataModelLib) $(ZMQLib) $(ZMQInclude) $(MyToolsInclude)  $(BoostLib) $(BoostInclude)


clean: 
	@echo -e "\n*************** Cleaning up ****************"
	rm -f include/*.h
	rm -f lib/*.so
	rm -f main
	rm -f RemoteControl
	rm -f NodeDaemon
	rm -f UserTools/*/*.o
	rm -f DataModel/*.o

lib/libDataModel.so: DataModel/* lib/libLogging.so lib/libStore.so $(patsubst DataModel/%.cpp, DataModel/%.o, $(wildcard DataModel/*.cpp))
	@echo -e "\n*************** Making " $@ "****************"
	cp DataModel/*.h include/
	#g++ -g -fPIC -shared DataModel/*.cpp -I include -L lib -lStore  -lLogging  -o lib/libDataModel.so $(DataModelInclude) $(DataModelLib) $(ZMQLib) $(ZMQInclude)  $(BoostLib) $(BoostInclude)
	g++  -fPIC -shared DataModel/*.o -I include -L lib -lStore -lLogging -o lib/libDataModel.so $(DataModelInclude) $(DataModelLib) $(ZMQLib) $(ZMQInclude) $(BoostLib) $(BoostInclude)

lib/libMyTools.so: UserTools/*/* UserTools/* include/Tool.h  lib/libLogging.so lib/libStore.so  $(patsubst UserTools/%.cpp, UserTools/%.o, $(wildcard UserTools/*/*.cpp)) |lib/libDataModel.so
	@echo -e "\n*************** Making " $@ "****************"
	cp UserTools/*/*.h include/
	cp UserTools/*.h include/
	#g++ -g -fPIC -shared  UserTools/Factory/Factory.cpp -I include -L lib -lStore -lDataModel -lLogging -o lib/libMyTools.so $(MyToolsInclude) $(MyToolsLib) $(DataModelInclude) $(DataModelLib) $(ZMQLib) $(ZMQInclude) $(BoostLib) $(BoostInclude)
	g++   -shared -fPIC UserTools/*/*.o -I include -L lib -lStore -lDataModel -lLogging -o lib/libMyTools.so $(MyToolsInclude) $(DataModelInclude) $(MyToolsLib) $(ZMQLib) $(ZMQInclude) $(BoostLib) $(BoostInclude)

RemoteControl:
	cd $(ToolDAQPath)/ToolDAQFramework/ && make RemoteControl
	@echo -e "\n*************** Copying " $@ "****************"
	cp $(ToolDAQPath)/ToolDAQFramework/RemoteControl ./

NodeDaemon: 
	cd $(ToolDAQPath)/ToolDAQFramework/ && make NodeDaemon
	@echo -e "\n*************** Copying " $@ "****************"
	cp $(ToolDAQPath)/ToolDAQFramework/NodeDaemon ./

lib/libServiceDiscovery.so: $(ToolDAQPath)/ToolDAQFramework/src/ServiceDiscovery/* | lib/libStore.so
	cd $(ToolDAQPath)/ToolDAQFramework && make lib/libServiceDiscovery.so
	@echo -e "\n*************** Copying " $@ "****************"
	cp $(ToolDAQPath)/ToolDAQFramework/src/ServiceDiscovery/ServiceDiscovery.h include/
	cp $(ToolDAQPath)/ToolDAQFramework/lib/libServiceDiscovery.so lib/
	#g++ -shared -fPIC -I include $(ToolDAQPath)/ToolDAQFramework/src/ServiceDiscovery/ServiceDiscovery.cpp -o lib/libServiceDiscovery.so -L lib/ -lStore  $(ZMQInclude) $(ZMQLib) $(BoostLib) $(BoostInclude)

lib/libLogging.so:  $(ToolDAQPath)/ToolDAQFramework/src/Logging/* | lib/libStore.so
	cd $(ToolDAQPath)/ToolDAQFramework && make lib/libLogging.so
	@echo -e "\n*************** Copying " $@ "****************"
	cp $(ToolDAQPath)/ToolDAQFramework/src/Logging/Logging.h include/
	cp $(ToolDAQPath)/ToolDAQFramework/lib/libLogging.so lib/
	#g++ -shared -fPIC -I include $(ToolDAQPath)/ToolDAQFramework/src/Logging/Logging.cpp -o lib/libLogging.so -L lib/ -lStore $(ZMQInclude) $(ZMQLib) $(BoostLib) $(BoostInclude)

update:
	@echo -e "\n*************** Updating ****************"
	cd $(ToolDAQPath)/ToolDAQFramework; git pull
	cd $(ToolDAQPath)/zeromq-4.0.7; git pull
	git pull


UserTools/%.o: UserTools/%.cpp lib/libStore.so include/Tool.h lib/libLogging.so lib/libDataModel.so
	@echo -e "\n*************** Making " $@ "****************"
	cp $(shell dirname $<)/*.h include
	-g++ $(CXXFLAGS) -c -fPIC -o $@ $< -I include -L lib -lStore -lDataModel -lLogging $(MyToolsInclude) $(MyToolsLib) $(DataModelInclude) $(DataModelLib) $(ZMQLib) $(ZMQInclude) $(BoostLib) $(BoostInclude)

target: remove $(patsubst %.cpp, %.o, $(wildcard UserTools/$(TOOL)/*.cpp))

remove:
	echo "removing tool objects"
	-rm UserTools/$(TOOL)/*.o

DataModel/%.o: DataModel/%.cpp lib/libLogging.so lib/libStore.so
	@echo -e "\n*************** Making " $@ "****************"
	cp $(shell dirname $<)/*.h include
	-g++ $(CXXFLAGS) -c -fPIC -o $@ $< -I include -L lib -lStore -lLogging  $(DataModelInclude) $(DataModelLib) $(ZMQLib) $(ZMQInclude) $(BoostLib) $(BoostInclude)

# this is the 'hack' library that must be linked into the final executable,
# which together with certain orderings of arguments (place the library EARLY)
# somehow avoids cernlib 64-bit address issues
lib/liblowfit_sk4_stripped.so: DataModel/lowfit_sk4_stripped.cc DataModel/lowfit_sk4_stripped.h
	@echo -e "\n*************** Making " $@ "****************"
	-g++ $(CXXFLAGS) -shared -L${CERN_ROOT}/lib $< -o $@

