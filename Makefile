# 设置编译器
FC = ifort

# 设置源文件
SRCS = n3lo500new.f deuteron.f90

# 设置输出文件名
OUT = deuteron

# 编译规则
all: $(OUT)

$(OUT): $(SRCS)
	$(FC) $(SRCS) -o $(OUT) -llapack

clean:
	rm -f $(OUT) *.mod