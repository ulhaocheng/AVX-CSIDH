# AVXCSIDH

AVXCSIDH is a software library that includes the batched high-throughput and the
unbatched low-latency AVX-512 implementations of CSIDH.

## The Structure of AVX-CSIDH Software Library
<pre> AVX-CSIDH
      ├── AVX-512IFMA implementations 
      |      ├── OAYT-style 
      |      |     ├── unbatched low-latency implementation 
      |      |     └── batched high-throughput implementation
      |      |          ├── extra-dummy    batching method 
      |      |          ├── extra-infinity batching method
      |      |          └── combined       batching method 
      |      └── DummyFree-style 
      |            ├── unbatched low-latency implementation 
      |            └── batched high-throughput implementation
      |                 ├── extra-dummy    batching method 
      |                 ├── extra-infinity batching method 
      |                 └── combined       batching method 
      |
      └── AVX-512F implementations 
             ├── OAYT-style 
             |     └── batched high-throughput implementation
             |          ├── extra-dummy    batching method 
             |          ├── extra-infinity batching method
             |          └── combined       batching method 
             └── DummyFree-style 
                   └── batched high-throughput implementation
                        ├── extra-dummy    batching method 
                        ├── extra-infinity batching method 
                        └── combined       batching method 
</pre>

## Usage

### Compile the AVX-512IFMA implementation (need Intel Cannon Lake, Ice Lake or Tiger Lake machine!)

**The batched high-throughput (ht) implementation:** 
```bash
    $ cd AVX-CSIDH/AVX-512IFMA-version 
    $ make ht_csidh STYLE=[OAYT/DUMMY_FREE] METHOD=[EXTRA_DUMMY/EXTRA_INFINITY/COMBINED]
    $ ./ht_csidh
```
For example, if you'd like to compile the AVX-512IFMA OAYT-style high-throughput
implementation which uses the combined batching method, then the 2nd line is:
```make ht_csidh STYLE=OAYT METHOD=COMBINED```

**The unbatched low-latency (ll) implementation:** 
```bash
    $ cd AVX-CSIDH/AVX-512IFMA-version 
    $ make ll_csidh STYLE=[OAYT/DUMMY_FREE] 
    $ ./ll_csidh
```

For example, if you'd like to compile the AVX-512IFMA DummyFree-style low-latency implementation, then the 2nd line is: 
```make ll_csidh STYLE=DUMMY_FREE``` 

### Compile the AVX-512F implementation (need Intel AVX-512 machine!)

```bash
    $ cd AVX-CSIDH/AVX-512F-version 
    $ make ht_csidh STYLE=[OAYT/DUMMY_FREE] METHOD=[EXTRA_DUMMY/EXTRA_INFINITY/COMBINED]
    $ ./ht_csidh
```
For example, if you'd like to compile the AVX-512F OAYT-style high-throughput
implementation which uses the extra-dummy batching method, then the 2nd line is: 
```make ht_csidh STYLE=OAYT METHOD=EXTRA_DUMMY```

## Paper
An paper describing the various implementations in this library has been
published in *IACR Transactions on Cryptographic Hardware and Embedded Systems,
2021(4), 618-649*.

The paper online link is
[here](https://tches.iacr.org/index.php/TCHES/article/view/9077).

## References 
  * [CCC+19] D. Cervantes-Vázquez, M. Chenu, J.-J. Chi-Domínguez, L. De Feo, F.
    Rodríguez Henríquez, and B. Smith. *Stronger and faster side-channel
    protections for csidh.* In P. Schwabe and N. Thériault, editors, Progress in
    Cryptology – LATINCRYPT 2019, pages 173–193. Springer International
    Publishing, 2019.

## Software Author
Hao Cheng (University of Luxembourg).

## Copyright
Copyright (C) 2021 by University of Luxembourg.

## LICENSE
GPLv3 (see details in LICENSE file).
