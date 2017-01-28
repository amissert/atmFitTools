{

gROOT->ProcessLine(".L markovDiff.cxx+");

//markovDiff* md = new markovDiff("mcmctree.root",1);
//markovDiff* md = new markovDiff("./mcmctree.root",10000);
markovDiff* md = new markovDiff("./run/results/fit*.root",10000,1);

}
