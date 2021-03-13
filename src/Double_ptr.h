class Double_ptr
{
public:
    Double_ptr(size_t n);
    ~Double_ptr();
    double* get_ptr();
    double* get_another_ptr();
    void switch_ptr();
    size_t length;
private:
    bool is_switched;
    double *ptr1;
    double *ptr2;
};

