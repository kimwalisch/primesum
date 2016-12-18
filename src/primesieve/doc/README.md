Build manpage & API documentation
=================================

In order to build the primesieve manpage and the html API
documentation you need to have installed the ```help2man``` and
```doxygen``` programs. Run the commands below from the parent
directory.

```bash
./configure --enable-maintainer-mode
make man
make doxygen
```
