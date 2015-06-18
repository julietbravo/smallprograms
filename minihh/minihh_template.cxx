#include <iostream>

template <class dtype>
class Advec
{
    public:
        Advec();
        void exec(dtype*, int);
};

template <class dtype>
class Field
{
    public:
        Field(int);
        int size;
        dtype* data; 
};

template <class dtype>
class Model
{
    public:
        Model();
        void run();

    private:
        Advec<dtype>* advec;
        Field<dtype>* field;
};

template <class dtype>
Advec<dtype>::Advec()
{
    std::cout << "creating Advec" << std::endl;
}

template <class dtype>
Field<dtype>::Field(int n)
{
    std::cout << "creating Field" << std::endl;
    data = new dtype[n];
    size = n;
}

template <class dtype>
void Advec<dtype>::exec(dtype* data, int n)
{
    std::cout << "exec Advec" << std::endl;
    std::cout << "dtype size=" << sizeof(data[0]) << " bytes" << std::endl;
}

template <class dtype>
Model<dtype>::Model()
{
    std::cout << "creating Model" << std::endl;

    const int n = 20;

    field = new Field<dtype>(n);
    advec = new Advec<dtype>;
}

template <class dtype>
void Model<dtype>::run()
{
    advec->exec(field->data, field->size);
}

int main()
{
    Model<float> model;

    model.run();

    return 0;
}
