# ASAGI I/O

[ASAGI](https://github.com/TUM-I5/ASAGI) is a fileformat that is used by codes such as SeisSol or ExaHype. It employs the NetCDF4 data format.

We can read ASAGI files into GMG, and write GMG datasets to ASAGI format, which makes it straightforward to use GMG to create a model setup or to couple the results of geodynamic simulations (e.g., produced by LaMEM) with codes that support ASAGI.

```@docs
read_ASAGI
write_ASAGI
```
