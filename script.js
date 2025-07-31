const canvas = document.getElementById('starfield');
const ctx = canvas.getContext('2d');

function resizeCanvas() {
    canvas.width = window.innerWidth;
    canvas.height = window.innerHeight;
}

window.addEventListener('resize', resizeCanvas);
resizeCanvas();

const mouse = {
    x: canvas.width / 2,
    y: canvas.height / 2
};

window.addEventListener('mousemove', (event) => {
    mouse.x = event.clientX;
    mouse.y = event.clientY;
});

class Star {
    constructor() {
        this.x = Math.random() * canvas.width;
        this.y = Math.random() * canvas.height;
        this.z = Math.random() * canvas.width;
        this.size = Math.random() * 2 + 0.5;
        this.speedX = (Math.random() - 0.5) * 0.5;
        this.speedY = (Math.random() - 0.5) * 0.5;
    }

    update() {
        this.z -= 1;
        if (this.z < 1) {
            this.z = canvas.width;
            this.x = Math.random() * canvas.width;
            this.y = Math.random() * canvas.height;
        }

        const dx = this.x - mouse.x;
        const dy = this.y - mouse.y;
        const distance = Math.sqrt(dx * dx + dy * dy);
        const forceDirectionX = dx / distance;
        const forceDirectionY = dy / distance;
        const maxDistance = 200;
        const force = (maxDistance - distance) / maxDistance;

        if (distance < maxDistance) {
            this.x += forceDirectionX * force * 2;
            this.y += forceDirectionY * force * 2;
        }

        this.x += this.speedX;
        this.y += this.speedY;

        if (this.x > canvas.width) this.x = 0;
        if (this.x < 0) this.x = canvas.width;
        if (this.y > canvas.height) this.y = 0;
        if (this.y < 0) this.y = canvas.height;
    }

    draw() {
        const x = (this.x - canvas.width / 2) * (canvas.width / this.z) + canvas.width / 2;
        const y = (this.y - canvas.height / 2) * (canvas.width / this.z) + canvas.height / 2;
        const radius = this.size * (canvas.width / this.z);

        ctx.beginPath();
        ctx.arc(x, y, radius, 0, Math.PI * 2);
        ctx.fillStyle = 'white';
        ctx.fill();
    }
}

const stars = [];
for (let i = 0; i < 500; i++) {
    stars.push(new Star());
}

function animate() {
    ctx.clearRect(0, 0, canvas.width, canvas.height);
    stars.forEach(star => {
        star.update();
        star.draw();
    });
    requestAnimationFrame(animate);
}

animate();