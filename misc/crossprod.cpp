using Eigen::MatrixXd;
using Eigen::Upper;
typedef Eigen::Map<MatrixXd>  MMat;

MMat V(as<MMat> (RV));
MMat VtV(as<MMat> (RVtV));
printf("V.data(): %p\n", (void *) V.data());
printf("V dims: %ld, %ld\n", V.rows(), V.cols());
printf("VtV.data(): %p\n", (void *) VtV.data());
printf("VtV dims: %ld, %ld\n", VtV.rows(), VtV.cols());

VtV.setZero().selfadjointView<Upper>().rankUpdate(V.adjoint());

return wrap(VtV);
