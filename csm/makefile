# 定义编译器和链接器
FC = gfortran
LD = $(FC)
#定义一些git bash的shell命令变量
BASE := $(shell expr $(CURDIR) : "\(.*\)/.*")
VERDATE := $(shell git log -1 --format=%cd  )
VERREV := $(shell git log -1 --pretty=format:"%h")
COMPDATE :=$(shell date)
# 定义编译选项
FFLAGS = -O3 -Wall -Wextra  -ffixed-line-length-0
LDFLAGS = -llapack 
#LDFLAGS = -llapack -lgsl
# 定义需要链接的库
LIBS = -llapack 
#LIBS = -llapack -lgsl
# 定义目标文件
TARGET = csm
# 定义源文件
SOURCES = module.f variable.f parameter.f eigenvalue.f algorithm.f pot.f hobasis.f thobasis.f csm.f
# 将源文件转换为对象文件
OBJECTS = $(SOURCES:.f=.o)

# 默认目标
all: $(TARGET)

# 生成可执行文件
$(TARGET): $(OBJECTS)
	$(LD) $(LDFLAGS) $(OBJECTS) -o $(TARGET) $(LIBS)

# 生成对象文件
%.o: %.f %.f90 %.f95
	$(FC) $(FFLAGS) -c $< -o $@

# 清理中间文件和可执行文件
clean:
	rm -f $(OBJECTS) $(TARGET) *.mod fort.*

.PHONY: print-info

print-info:
	@echo "Base path: $(BASE)"
	@echo "Version Date: $(VERDATE)"
	@echo "Version Revision: $(VERREV)"
	@echo "Compilation Date: $(COMPDATE)"
