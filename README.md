# Zero-Knowledge Proof Learning Resources (a learning path)

Zero-knowledge proof (ZKP) probably is currently the hottest topic in the cryptography world. The idea of proving that `something is true without revealing any information apart from the fact that this specific statement is true` is one of the most beautiful things that cryptography could do.


Although the [ZKP concept](http://groups.csail.mit.edu/cis/pubs/shafi/1985-stoc.pdf) was introduced a long time ago by MIT researchers Shafi Goldwasser, Silvio Micali, and Charles Rackoff, it was just getting more attention to practical applications recently and yields a huge potential in Web3. 

From my experience, I found that it is hard to understand zero-knowledge through academic papers. If you want, I think the most comprehensive explanation trying to illustrate this concept was probably written by Matthew Green in his blog and and YouTube video by Avi Wigderson:

- [Zero Knowledge Proofs: An illustrated primer](https://blog.cryptographyengineering.com/2014/11/27/zero-knowledge-proofs-illustrated-primer/), written by Matthew Green.
- [Zero Knowledge Proof - Numberphile](https://www.youtube.com/watch?v=5ovdoxnfFVc&ab_channel=Numberphile2), given by Avi Wigderson.

After getting a sense of what is zero-knowledge proof, you then can go to very good video lectures that summarize the development of zero-knowledge proof schemes:

- [Introduction to Zero Knowledge Interactive Proofs](https://youtube.com/playlist?list=PLS01nW3RtgorR09s4cIz3aFylYCrk8fv0&si=77NRiOm7pkLk4OeX), given by Shafi Goldwasser. 
- [Overview of Modern SNARK Constructions](https://www.youtube.com/watch?v=bGEXYpt3sj0&ab_channel=Blockchain-Web3MOOCs), given by Dan Boneh. 
- [Bilinear Pairings-based Zero-Knowledge Proofs](https://www.youtube.com/watch?v=X-z3JYlFdzs&ab_channel=ZKProofStandards), given by Jens Groth.

That is probably enough to get started with ZKP from theory. Recently, when you heard about Zero-Knowledge systems on BLockchain or Web3, that is more about `applied` ZKP. Thanks to some ZK languages such as [circom](https://docs.circom.io/), [halo2](https://zcash.github.io/halo2/), [noir](https://aztec.network/noir/), etc., your task to implement a ZP program is much easier. 

- [The RareSkills Book of Zero Knowledge](https://www.rareskills.io/zk-book)
- [Zero-Knowledge Proofs MOOC](https://www.youtube.com/@blockchain-web3moocs635#)

Circom and Halo2 are probably the most popular ZK languages. Their documentation is in the following links:

- [Step-by-step: how to implement and deploy a Circom's circuit](https://docs.circom.io/)
- [The Halo2 book](https://zcash.github.io/halo2), by ZCash
- [0xPARC](https://learn.0xparc.org/): a fruitful resources introducing the use of Circom and halo2.


Personally, I found that the Rareskills book is great source. It explained and demonstrated by Python step-by-step how to a practical Zero-Knowledge Proof system work. You can find code in the subfolder `Python`. 
