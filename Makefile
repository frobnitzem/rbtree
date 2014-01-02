OBJS=rbtree.o

rbtree.a:	$(OBJS)
	$(AR) -cr $@ $^

test:	rbtree.a test.o
	$(LD) -o $@ $^

distclean: clean
	rm -f rbtree.a

clean:
	rm -f $(OBJS)

.c.o:	$(PWD)/include
	$(CC) $(CFLAGS) -c -o $@ $^

