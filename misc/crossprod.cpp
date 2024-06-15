using Eigen::MatrixXd;
using Eigen::Upper;
typedef Eigen::Map<MatrixXd>  MMat;

MMat A(as<MMat> (AA));
const int m(A.rows()), n(A.cols());
MatrixXd firstMat(n,n);
printf("firstMat.data(): %p\n", (void *) firstMat.data());
printf("firstMat dims: %ld, %ld\n", firstMat.rows(), firstMat.cols());
MatrixXd AtA(firstMat.setZero().selfadjointView<Upper>().rankUpdate(A.adjoint()));
return wrap(AtA);
