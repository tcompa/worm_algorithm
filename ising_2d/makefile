CC:=g++ -g -Wall -lm

worm_ising_2d: worm_ising_2d.o
	${CC} -o $@.x $^  -lm

%.o : %.cpp
	${CC} -c $<

clean:
	rm *.o *.x *~
