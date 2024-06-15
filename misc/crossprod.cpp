using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::Lower;
using Eigen::Upper;

const Map<MatrixXd> A(as<Map<MatrixXd> >(AA));
const int m(A.rows()), n(A.cols());
MatrixXd firstMat(n,n);
printf("firstMat.data(): %p\n", (void *) firstMat.data());
printf("firstMat dims: %ld, %ld\n", firstMat.rows(), firstMat.cols());
MatrixXd AtA(firstMat.setZero().selfadjointView<Upper>().rankUpdate(A.adjoint()));
MatrixXd AAt(MatrixXd(m, m).setZero().selfadjointView<Upper>().rankUpdate(A));
return List::create(Named("crossprod(A)") = AtA,
		    Named("tcrossprod(A)") = AAt);
