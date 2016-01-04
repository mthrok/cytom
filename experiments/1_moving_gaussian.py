import numpy as np
import sklearn.mixture


def parse_command_line_arguments():
    from argparse import ArgumentParser as AP
    ap = AP()
    ap.add_argument('--sigma', type=float, default=1.0)
    ap.add_argument('--n-samples', type=int, default=1000)
    ap.add_argument('--n-clusters', type=int, default=10)
    ap.add_argument('--vx', type=float, default=3.0)
    ap.add_argument('--vy', type=float, default=0.0)
    return ap.parse_args()


def generate_mus(vx, vy, n_steps=10):
    mus = np.zeros((n_steps, 2))
    for i in range(n_steps - 1):
        mus[i + 1] = mus[i] + [vx, vy]
    return mus


def generate_data(sigma, n_samples, vx, vy, n_steps=10):
    mus = generate_mus(vx, vy, n_steps=n_steps)
    cov = sigma * np.identity(2)
    data = np.zeros((n_samples * n_steps, 2))
    for i, mu in enumerate(mus):
        d = np.random.multivariate_normal(mu, cov, n_samples)
        data[i*n_samples:(i+1)*n_samples] = d
    return {'mus': mus, 'data': data}


def plot_data(mus, data, labels, n_clusters):
    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax = fig.add_subplot(2, 2, 2)
    for l in range(n_clusters):
        index = (labels == l)
        color = np.random.rand(3, 1)
        ax.scatter(data[index, 0], data[index, 1], color=color, s=1)
        ax.scatter(np.mean(data[index, 0]), np.mean(data[index, 1]))
    ax.scatter(mus[:, 0], mus[:, 1], color='red')
    ax.set_aspect('equal', adjustable='box')
    ax = fig.add_subplot(2, 2, 3)
    ax.hexbin(data[:, 0], data[:, 1], gridsize=20)
    ax.set_aspect('equal', adjustable='box')
    ax = fig.add_subplot(2, 2, 1)
    ax.hist(data[:, 0])
    ax = fig.add_subplot(2, 2, 4)
    ax.hist(data[:, 1], orientation='horizontal')
    plt.show()


def analyze_data(data, n_components):
    gmm = sklearn.mixture.GMM(n_components=n_components)
    gmm.fit(data)
    labels = gmm.predict(data)
    print(gmm.weights_)
    print(gmm.means_)
    print(gmm.covars_)
    return labels


def main():
    args = parse_command_line_arguments()
    data = generate_data(args.sigma, args.n_samples, args.vx, args.vy)
    labels = analyze_data(data['data'], args.n_clusters)
    plot_data(data['mus'], data['data'], labels, args.n_clusters)

if __name__ == '__main__':
    main()
