template<class T>
class BinaryTree {
public:
    typedef typename T::Face Face;
    typedef std::vector<Face> Faces;

    class Node {
    public:
        Faces data;
        Node *left;
        Node *right;

        // Val is the key or the value that
        // has to be added to the data part
        explicit Node(Faces val) {
            data = val;

            // Left and right child for node
            // will be initialized to null
            left = nullptr;
            right = nullptr;
        }

        void printData() {
            std::cout << "[";
            for (auto var: data) {
                typename T::Triangle a(var.vertex(0)->point(), var.vertex(1)->point(), var.vertex(2)->point());
                std::cout << a;
            }
            std::cout << "]" << std::endl;
        }
    };


    Node *Root() {
        return root;
    }

    void printPreOrder(Node *node) {
        if (node == nullptr)
            return;

        for (auto var: node->data) {
            std::cout << var.vertex(0)->point() << var.vertex(1)->point() << var.vertex(2)->point() << std::endl;
        }
        printPreOrder(node->left);
        printPreOrder(node->right);
    }

    void printInOrder(Node *node) {
        if (node == nullptr)
            return;

        printInOrder(node->left);
        for (auto var: node->data) {
            std::cout << var.vertex(0)->point() << var.vertex(1)->point() << var.vertex(2)->point() << std::endl;
        }
        printInOrder(node->right);
    }

    void printPostOrder(Node *node) {
        if (node == nullptr)
            return;

        printPostOrder(node->left);
        printPostOrder(node->right);
        node->printData();
    }

    void prettyPrint() {
        if (root) {
            printHelper(this->root, "", true);
        }
    }

    BinaryTree(Faces val) {
        root = new Node(val);
    }

    BinaryTree() {};

    Node *EmptyNode() {
        return new Node(Faces());
    }

private:
    Node *root;

    void printHelper(Node *rootArg, std::string indent, bool last) {
        // print the tree structure on the screen
        if (rootArg != nullptr) {
            std::cout << indent;
            if (last) {
                std::cout << "R----";
                indent += "     ";
            } else {
                std::cout << "L----";
                indent += "|    ";
            }

//            std::cout << rootArg->data.size() << std::endl;
            rootArg->printData();
            printHelper(rootArg->left, indent, false);
            printHelper(rootArg->right, indent, true);
        }
    }
};
